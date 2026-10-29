"""Apply the released pertpy Mixscape to fixed, verified single-guide cohorts.

Output labels describe RNA response under the model, never DNA editing or
functional knockout efficiency. All design guides remain in the coverage table.
The official signature and classifier are called without reimplementation.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd

from compare_established_perturbseq import load_counts, sha_file

SCIENTIFIC_CONTRACT = dict(package='pertpy', version='1.4.0', total_RNA_scale=10000,
                           HVG_flavor='seurat', HVGs=2000, PCA=40, NTC_neighbors=20,
                           min_DE_RNAs=5, logfc_threshold=0.25, adjusted_p_cutoff=0.05,
                           iterations=10, seed=20261029)


def atomic(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf8')
    temporary.replace(path)


def model_labels(cells):
    """Pool only genuine target identities; keep intergenic guides separate."""
    result = cells.target.astype(str).copy()
    nt = cells.role.eq('non_targeting_control')
    intergenic = cells.role.eq('intergenic_control')
    result.loc[nt] = 'NT'
    result.loc[intergenic] = 'INTERGENIC:' + cells.loc[intergenic, 'guide_id'].astype(str)
    if result.eq('').any() or result.isna().any():
        raise ValueError('Missing model group identity')
    return result


def summarize_guides(catalogue, cells, classified):
    """Missing measurements stay missing, including the controls themselves."""
    if not catalogue.guide.is_unique or not cells.cell_barcode.is_unique:
        raise ValueError('Duplicate guide or cell identities')
    if not classified.index.is_unique or not classified.index.isin(cells.cell_barcode).all():
        raise ValueError('Classification axes do not match input cells')
    joined = cells.set_index('cell_barcode').join(classified, how='left')
    output = catalogue[['guide', 'target_gene_name', 'source_control_type']].copy()
    counts = joined.groupby('guide_id').size()
    evaluated = joined.rna_class.isin(['KO', 'NP'])
    n_evaluated = joined[evaluated].groupby('guide_id').size()
    n_response = joined[joined.rna_class.eq('KO')].groupby('guide_id').size()
    output['n_QC_records'] = output.guide.map(counts).fillna(0).astype(int)
    output['n_model_evaluated'] = output.guide.map(n_evaluated).fillna(0).astype(int)
    output['n_model_RNA_response'] = output.guide.map(n_response).fillna(0).astype(int)
    output['model_RNA_response_fraction'] = np.divide(
        output.n_model_RNA_response, output.n_model_evaluated,
        out=np.full(len(output), np.nan), where=output.n_model_evaluated.to_numpy() > 0)
    output['status'] = 'evaluated_RNA_model_not_DNA_KO'
    output.loc[output.n_model_evaluated.eq(0), 'status'] = 'insufficient_target_group_cells'
    output.loc[output.source_control_type.eq('non_targeting_control'), 'status'] = 'NTC_reference_not_scored'
    output.loc[output.n_QC_records.eq(0), 'status'] = 'no_QC_records'
    return output, joined


def apply_official(a):
    """Call both stages of released Mixscape with the frozen public settings."""
    import pertpy as pt
    ms = pt.tl.Mixscape()
    ms.perturbation_signature(a, 'model_group', 'NT', n_neighbors=20,
                             use_rep='X_pca', n_dims=40, split_by=None,
                             ref_selection_mode='nn', random_state=20261029, n_jobs=2)
    ms.mixscape(a, 'model_group', 'NT', layer='X_pert', min_de_genes=5,
                logfc_threshold=0.25, pval_cutoff=0.05, test_method='wilcoxon',
                iter_num=10, scale=True, split_by=None, random_state=20261029)


def checkpoint_strings(a):
    """Use supported string storage without changing labels or row order."""
    for frame in [a.obs, a.var]:
        frame.index = pd.Index(frame.index.to_numpy(dtype=object), dtype=object)
        for column in frame:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = pd.Series(frame[column].to_numpy(dtype=object),
                                          index=frame.index, dtype=object)
    return a


def ntc_holdout_diagnostic(a, destination):
    """Evaluate whole NTC guides excluded from the reference, not cell halves.

    The expression representation is frozen. This descriptive negative-control
    check is not a formal false-positive calibration or a DNA genotype test.
    Its results never select controls or change the main method's thresholds.
    """
    controls = a[a.obs.role.eq('non_targeting_control')].copy()
    guides = sorted(controls.obs.guide_id.unique(),
                    key=lambda g: hashlib.sha256(str(g).encode()).hexdigest())
    query = guides[::5]
    held = controls.obs.guide_id.isin(query)
    controls.obs['model_group'] = 'NT'
    controls.obs.loc[held, 'model_group'] = 'NTC_QUERY:' + controls.obs.loc[held, 'guide_id'].astype(str)
    count = controls.obs.model_group.value_counts()
    controls = controls[controls.obs.model_group.map(count).ge(2)].copy()
    if controls.obs.model_group.eq('NT').sum() <= 20:
        raise ValueError('Insufficient disjoint controls for held-guide diagnostic')
    apply_official(controls)
    controls.obs.to_csv(destination / 'NTC_held_guide_classifications.csv.gz')
    rows = []
    for guide in query:
        q = controls.obs[controls.obs.guide_id.eq(guide)]
        n = len(q)
        k = int(q.mixscape_class_global.eq('KO').sum()) if n else 0
        rows.append(dict(guide=guide, n_evaluated=n, n_model_response=k,
                         model_response_fraction=k/n if n else None))
    pd.DataFrame(rows).to_csv(destination / 'NTC_held_guide_summary.csv', index=False)
    atomic(destination / 'NTC_reference_isolation.json', dict(
        queried_whole_guides=query,
        reference_whole_guides=[g for g in guides if g not in query],
        selection='every fifth guide after SHA256(name) ordering; no RNA-effect selection',
        query_reference_guide_overlap=0, reference_unchanged_for_main_analysis=True,
        formal_false_positive_calibration=False))
    return dict(queried_guides=len(query), guides_with_any_model_response=sum(r['n_model_response'] > 0 for r in rows),
                query_records=sum(r['n_evaluated'] for r in rows),
                query_records_model_response=sum(r['n_model_response'] for r in rows))


def run_dataset(spec, destination, progress):
    import anndata as ad
    import scanpy as sc
    import pertpy as pt

    destination.mkdir(parents=True, exist_ok=True)
    checkpoint = destination / 'normalized_HVG_checkpoint.h5ad'
    receipt = destination / 'checkpoint_receipt.json'
    config_sha = hashlib.sha256(json.dumps(dict(dataset=spec, method=SCIENTIFIC_CONTRACT),
                                          sort_keys=True).encode()).hexdigest()
    completed = destination / 'summary.json'
    manifest = destination / 'file_manifest.json'
    if completed.exists() and manifest.exists():
        result = json.loads(completed.read_text())
        if result['spec_sha256'] != config_sha:
            raise ValueError('Completed dataset contract differs')
        saved_files = json.loads(manifest.read_text())
        for name, item in saved_files.items():
            file = destination / name
            if not file.exists() or file.stat().st_size != item['bytes']:
                raise ValueError('Completed result missing or changed')
            if name != checkpoint.name and sha_file(file) != item['sha256']:
                raise ValueError('Completed small result hash differs')
        progress('reuse_completed_dataset', dataset=spec['id'])
        return result
    if checkpoint.exists() or receipt.exists():
        if not (checkpoint.exists() and receipt.exists()):
            raise ValueError('Incomplete checkpoint publication; preserve for diagnosis')
        saved = json.loads(receipt.read_text())
        if saved['config_sha256'] != config_sha or sha_file(checkpoint) != saved['sha256']:
            raise ValueError('Checkpoint identity or hash differs')
        progress('resume_verified_HVG_checkpoint', dataset=spec['id'])
        a = ad.read_h5ad(checkpoint)
        cells = pd.read_csv(destination / 'all_input_cells.csv.gz', dtype={'cell_barcode': str})
        catalogue = pd.read_csv(spec['input']['catalogue'], keep_default_na=False)
        if sha_file(destination / 'all_input_cells.csv.gz') != saved['cells_sha256']:
            raise ValueError('Checkpoint cell metadata differs')
    else:
        progress('read_verified_counts', dataset=spec['id'])
        x, cells, genes, proof = load_counts(spec['input'], progress)
        catalogue = pd.read_csv(spec['input']['catalogue'], keep_default_na=False)
        cells['model_group'] = model_labels(cells)
        group_n = cells.model_group.value_counts()
        # Scanpy's test requires at least two observations in each tested group.
        keep = cells.model_group.map(group_n).ge(2).to_numpy()
        if (cells.role.eq('non_targeting_control') & keep).sum() <= 20:
            raise ValueError('Insufficient NTCs for the official 20-neighbor signature')
        obs = cells.loc[keep].set_index('cell_barcode').copy()
        obs.index = obs.index.astype(str)
        a = ad.AnnData(x[keep].astype(np.float32), obs=obs,
                       var=genes.set_index('gene_id').copy())
        if not a.var_names.is_unique:
            raise ValueError('RNA feature IDs must be unique')
        sc.pp.normalize_total(a, target_sum=10000)
        sc.pp.log1p(a)
        # The released pertpy tutorial subsets to Seurat-flavor HVGs before
        # the official signature/classification; freeze 2000, not an r filter.
        sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor='seurat', subset=True)
        sc.pp.pca(a, n_comps=40, random_state=20261029)
        checkpoint_strings(a)
        cells.to_csv(destination / 'all_input_cells.csv.gz', index=False)
        atomic(destination / 'count_input_proof.json', proof)
        a.var.to_csv(destination / 'tested_RNA_features.csv')
        temporary = checkpoint.with_name(checkpoint.name + '.partial')
        a.write_h5ad(temporary, compression='gzip')
        temporary.replace(checkpoint)
        atomic(receipt, dict(config_sha256=config_sha, sha256=sha_file(checkpoint),
                            cells_sha256=sha_file(destination / 'all_input_cells.csv.gz'),
                            records=a.n_obs, RNAs=a.n_vars, input_proof=proof,
                            input_scope=spec['input']['input_scope']))
    progress('whole_NTC_guide_holdout_diagnostic', dataset=spec['id'], records=a.n_obs, RNAs=a.n_vars)
    negative_check = ntc_holdout_diagnostic(a, destination)
    progress('official_Mixscape_classification', dataset=spec['id'])
    apply_official(a)
    classified = a.obs[['mixscape_class_global', 'mixscape_class_p_ko']].rename(
        columns={'mixscape_class_global': 'rna_class',
                 'mixscape_class_p_ko': 'model_RNA_response_probability'})
    if not classified.rna_class.isin(['NT', 'NP', 'KO']).all():
        raise ValueError('Unexpected official model label')
    guide_table, joined = summarize_guides(catalogue, cells, classified)
    guide_table.to_csv(destination / 'all_design_guides_RNA_model.csv', index=False)
    joined.to_csv(destination / 'all_input_cell_classifications.csv.gz')
    # Preserve the official first-iteration projection, including NTC values,
    # for checking separation. No arbitrary representative example is selected.
    projection = []
    for target, categories in a.uns['mixscape'].items():
        for split, values in categories.items():
            values = values.copy()
            values.insert(0, 'cell_barcode', values.index.astype(str))
            values.insert(1, 'target', target)
            values.insert(2, 'split', split)
            projection.append(values.reset_index(drop=True))
    if projection:
        pd.concat(projection, ignore_index=True).to_csv(
            destination / 'official_projection_distributions.csv.gz', index=False)
    targeting = guide_table.source_control_type.isin(['targeting', 'positive_control'])
    summary = dict(dataset=spec['id'], spec_sha256=config_sha, scope=spec['input']['input_scope'],
                   design_guides=len(guide_table), observed_guides=int(guide_table.n_QC_records.gt(0).sum()),
                   input_QC_records=len(cells), tested_RNAs=a.n_vars,
                   targeting_guides_model_evaluated=int((targeting & guide_table.n_model_evaluated.gt(0)).sum()),
                   targeting_records_model_response=int(guide_table.loc[targeting, 'n_model_RNA_response'].sum()),
                   targeting_records_model_evaluated=int(guide_table.loc[targeting, 'n_model_evaluated'].sum()),
                   targets_with_official_projection=len(a.uns['mixscape']),
                   held_whole_NTC_guide_diagnostic=negative_check,
                   DNA_editing_efficiency_assessed=False, functional_KO_efficiency_assessed=False,
                   negative_control_false_positive_calibration_completed=False,
                   limitation='NTC reference cells are not an independent scored negative-control experiment; '
                   'NP includes no detectable RNA response, not proof of failed DNA editing; '
                   'company merged barcode/capture identity limits remain.')
    atomic(destination / 'summary.json', summary)
    files = [p for p in destination.iterdir() if p.is_file() and not p.name.endswith('.partial')]
    checkpoint_sha = json.loads(receipt.read_text())['sha256']
    atomic(destination / 'file_manifest.json', {
        p.name: dict(bytes=p.stat().st_size,
                     sha256=checkpoint_sha if p == checkpoint else sha_file(p))
        for p in files})
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    output = Path(config['output'])
    output.mkdir(parents=True, exist_ok=True)
    lock = output / 'active.lock'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, json.dumps(dict(pid=os.getpid(), uid=os.getuid(), config=args.config)).encode())
    os.close(fd)
    start = time.monotonic()
    def progress(stage, **values):
        value = dict(stage=stage, at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                     elapsed_seconds=time.monotonic()-start, pid=os.getpid(), **values)
        atomic(output / 'status.json', value)
        print(json.dumps(value), flush=True)
    try:
        import pertpy as pt
        if importlib.metadata.version('pertpy') != '1.4.0':
            raise ValueError('Unexpected official package version')
        sources = {str(Path(inspect.getsourcefile(c))): sha_file(inspect.getsourcefile(c))
                   for c in [pt.tl.Mixscape, pt.tl.Mixscape.__bases__[0]]}
        atomic(output / 'software.json', dict(code_commit=config['code_commit'],
               official_package='pertpy', version='1.4.0', source_sha256=sources,
               classification_method='released official class; signature plus classifier',
               parameters=dict(normalize_total=10000, HVGs=2000, PCA=40,
                               NTC_neighbors=20, min_DE_RNAs=5, logfc_threshold=0.25,
                               adjusted_p_cutoff=0.05, iterations=10, seed=20261029)))
        summaries = []
        for spec in config['datasets']:
            summaries.append(run_dataset(spec, output / spec['id'], progress))
        atomic(output / 'summary.json', summaries)
        progress('complete_pending_scientific_interpretation')
    except BaseException as exc:
        progress('failed_pending_recovery', error=repr(exc))
        raise
    finally:
        lock.unlink()


if __name__ == '__main__':
    main()
