"""Challenge the adopted lab method on a selected public reference cohort.

The scientific functions are imported unchanged from the original lab pipeline.
Only the raw CSC input adapter and explicitly labelled diagnostic backgrounds
are new. The selected public release cannot validate its unavailable upstream
guide calling, library dropout or independent biological replication.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata as md
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

from analyze_guide_capture_qc import choose_qc
from audit_public_perturbseq import read_field
from perturbseq_core_io import atomic_json, identity, revision, select_features, sha256
from perturbseq_core_effects import (effect_profiles, expression_effect, pair_effects,
                                    qc_flow, stable_rng, target_support)
import perturbseq_core_plots as original_plots
import run_core_analysis as original_run
from public_count_cache import (CountMeans, axis_hash, cache_counts, moments,
                                selected_matrix)


REFERENCE_COMMIT = '993b6fab4e7c82c95f489bdf397274992eb9d5c4'
CONTRACT = dict(raw_mean_gt=.5, feature_cohort='additional_RNA_QC_single_guide',
                normalization='median_ntc_total_raw_RNA', min_guide_cells=10,
                matched_cells=10, matched_repeats=8, seed=20261022, pseudocount=1.,
                reference='all_QC_NTC_cells_except_whole_own_guide',
                batches='pooled_as_in_original_method',
                agreement='all_cells_per_guide_then_target_median_of_pair_Pearson',
                direct_RNA_mask='unique_symbol_in_selected_and_all_source_RNA')


def verify_original_functions(repo):
    """Byte equality of complete original science/plot modules, not a claim of validity."""
    proof = {}
    for name in ['perturbseq_core_effects.py', 'perturbseq_core_plots.py',
                 'run_core_analysis.py', 'analyze_guide_capture_qc.py']:
        original = subprocess.check_output(['git', 'show', REFERENCE_COMMIT+':scripts/'+name], cwd=repo)
        current = (Path(repo)/'scripts'/name).read_bytes().replace(b'\r\n', b'\n')
        if original.replace(b'\r\n', b'\n') != current:
            raise ValueError('Original scientific module changed: '+name)
        import hashlib
        proof[name] = hashlib.sha256(original).hexdigest()
    # The metadata reader evolved; pin the actual unchanged feature function.
    import ast
    old_io = subprocess.check_output(['git', 'show', REFERENCE_COMMIT+':scripts/perturbseq_core_io.py'], cwd=repo).decode()
    new_io = (Path(repo)/'scripts/perturbseq_core_io.py').read_text(encoding='utf8')
    def function_ast(text):
        return ast.dump(next(x for x in ast.parse(text).body
                             if isinstance(x, ast.FunctionDef) and x.name == 'select_features'))
    if function_ast(old_io) != function_ast(new_io):
        raise ValueError('Original RNA feature rule changed')
    return proof


def adapt_metadata(h, catalogue_path):
    keys = ['_index', 'gRNA', 'gene_target', 'batch', 'n_gRNA',
            'n_gRNA_UMIs', 'n_UMI_counts', 'n_genes']
    obs = pd.DataFrame({key: read_field(h['obs'][key]) for key in keys})
    genes = np.asarray(read_field(h['var/gene_ids']), dtype=str)
    names = np.asarray(read_field(h['var/_index']), dtype=str)
    if not obs['_index'].is_unique or not pd.Index(genes).is_unique or not obs.n_gRNA.eq(1).all():
        raise ValueError('Ambiguous cell/RNA identity or non-single-guide release')
    library = pd.read_csv(catalogue_path, dtype=str, keep_default_na=False)
    library['guide'] = library['id'].map(lambda value: '_'.join(value.rsplit('-', 1)))
    if library.guide.duplicated().any():
        raise ValueError('Ambiguous author guide-ID conversion')
    library = library.loc[library.guide.isin(obs.gRNA)].copy()
    mapping = library.set_index('guide').target_gene_name
    expected = obs.gene_target.replace({'NTC': 'Non-Targeting'}).astype(str)
    if len(library) != obs.gRNA.nunique() or not np.array_equal(obs.gRNA.map(mapping), expected):
        raise ValueError('Author guide/target catalogue differs from released cell assignments')
    library['source_control_type'] = np.where(library.target_gene_name.eq('Non-Targeting'),
                                            'non_targeting_control', 'targeting')
    obs.index = pd.Index(obs['_index'].astype(str), name='cell_barcode')
    obs['source_row'] = np.arange(len(obs), dtype=np.int64)
    obs['guide_id'] = obs.gRNA.astype(str)
    obs['guide_status'] = 'single_guide'
    obs['guide_umis'] = pd.to_numeric(obs.n_gRNA_UMIs, errors='raise')
    obs['raw_guide_umis'] = np.nan  # Not supplied. Never equate assigned and captured counts.
    var = pd.DataFrame(dict(gene_name=names), index=pd.Index(genes, name='gene_id'))
    return obs, var, library.sort_values('guide').reset_index(drop=True)


def diagnostic_backgrounds(obs, records, effects, pairs, seed):
    """One prespecified outcome-blind covariate match per real guide pair.

    These are descriptive comparators, not exchangeable independent null draws.
    No expression/correlation value is used to choose matches.
    """
    by_batch = pd.crosstab(obs.guide_id, obs.batch).reindex(records.guide, fill_value=0)
    fractions = by_batch.to_numpy(dtype=float, copy=True)
    fractions /= fractions.sum(1)[:, None]
    n = records.cells.to_numpy(float)
    guide_names = records.guide.to_numpy(str)
    target = records.target.to_numpy(str)
    roles = records.source_control_type.to_numpy(str)
    lookup = {g: i for i, g in enumerate(guide_names)}
    rows = []
    for pair in pairs.itertuples():
        a, b = lookup[pair.guide_a], lookup[pair.guide_b]
        cov_b = np.max(np.abs(fractions-fractions[b]), axis=1)
        count_b = np.abs(np.log(n/n[b]))
        candidates = np.flatnonzero((roles == 'targeting') & (target != target[a]) &
                                    (cov_b <= .1) & (count_b <= np.log(1.5)))
        if len(candidates):
            score = count_b[candidates] + np.abs(fractions[candidates]-fractions[b]).sum(1)
            order = np.lexsort((guide_names[candidates], score))
            # Deterministic nearest covariates; never choose a favorable r.
            c = int(candidates[order[0]])
            value = pair_effects(effects[a], effects[c])
            rows.append(dict(pair_target=pair.target, guide_a=pair.guide_a,
                             guide_b=pair.guide_b, comparator='different_target',
                             control_a=guide_names[a], control_b=guide_names[c],
                             matched=True, cells_a=int(n[a]), cells_b=int(n[c]),
                             max_batch_fraction_difference=float(cov_b[c]),
                             cell_count_ratio=float(n[c]/n[b]), **value))
        else:
            rows.append(dict(pair_target=pair.target, guide_a=pair.guide_a,
                             guide_b=pair.guide_b, comparator='different_target', matched=False))
        # Same tolerances for both sides of an NTC comparator.
        cov_a = np.max(np.abs(fractions-fractions[a]), axis=1)
        can_a = np.flatnonzero((roles == 'non_targeting_control') & (cov_a <= .1) &
                              (np.abs(np.log(n/n[a])) <= np.log(1.5)))
        can_b = np.flatnonzero((roles == 'non_targeting_control') & (cov_b <= .1) &
                              (count_b <= np.log(1.5)))
        options = [(int(x), int(y)) for x in can_a for y in can_b if x != y]
        if options:
            # Stable seed resolves equal metadata matches without effect selection.
            rng = stable_rng(seed, 'null:'+pair.guide_a+':'+pair.guide_b)
            x, y = options[int(rng.integers(len(options)))]
            value = pair_effects(effects[x], effects[y])
            rows.append(dict(pair_target=pair.target, guide_a=pair.guide_a,
                             guide_b=pair.guide_b, comparator='NTC_pair',
                             control_a=guide_names[x], control_b=guide_names[y], matched=True,
                             cells_a=int(n[x]), cells_b=int(n[y]),
                             max_batch_fraction_difference=float(max(cov_a[x], cov_b[y])),
                             cell_count_ratio=float(n[y]/n[b]), **value))
        else:
            rows.append(dict(pair_target=pair.target, guide_a=pair.guide_a,
                             guide_b=pair.guide_b, comparator='NTC_pair', matched=False))
    coverage = by_batch.copy()
    coverage.insert(0, 'guide', records.guide.to_numpy())
    coverage.index = np.arange(len(coverage))
    return pd.DataFrame(rows), coverage


def figure_annotations():
    """Presentation-only annotations for genuinely unavailable public inputs."""
    real_save = original_plots.save
    original_plots.NAMES['targeting'] = 'Targeting'
    original_run.NAMES['targeting'] = 'Targeting'

    def save(fig, path):
        name = Path(path).name
        if name == '01_technical_qc.png':
            fig.suptitle('Selected Strong release: additional RNA QC only; original guide-calling losses unavailable', fontsize=10)
        elif name == '02_four_type_coverage.png':
            fig.axes[0].set_title('Coverage of released sgRNAs after additional QC')
            fig.axes[2].clear()
            fig.axes[2].axis('off')
            fig.axes[2].text(.5, .5, 'Whole-library dropout: unavailable\nStrong is a selected release.\nIntergenic/positive roles not annotated.', ha='center', va='center', transform=fig.axes[2].transAxes)
        elif name == '03_four_type_cell_quality.png':
            fig.axes[3].clear(); fig.axes[3].axis('off')
            fig.axes[3].text(.5, .5, 'Captured guide UMIs\nnot supplied\nAssigned UMIs are shown separately', ha='center', va='center', transform=fig.axes[3].transAxes)
            fig.suptitle('Within released single-guide records; missing roles are unavailable', fontsize=10)
        real_save(fig, path)

    original_plots.save = save
    # This function imported save directly from the same unchanged original file.
    original_run.save = save


def run(config_path):
    config = json.loads(Path(config_path).read_text(encoding='utf8'))
    out = Path(config['output'])
    repo = Path(__file__).resolve().parents[1]
    code = revision(config['expected_commit'])
    proof = verify_original_functions(repo)
    source = Path(config['source'])
    if (not out.is_absolute() or out.parent != Path(config['allowed_output_parent']) or
            out.is_relative_to(repo) or source.is_relative_to(out)):
        raise ValueError('Distinct external public result location required')
    before = identity(source)
    if (before['bytes'], before['mtime_ns']) != (config['expected_bytes'], config['expected_mtime_ns']):
        raise ValueError('Accepted source identity changed')
    if sha256(config['catalogue']) != config['catalogue_sha256']:
        raise ValueError('Pinned author catalogue changed')
    out.mkdir(exist_ok=False)
    for folder in ['tables', 'figures', 'raw_count_cache']:
        (out/folder).mkdir()
    owner = dict(pid=os.getpid(), uid=os.getuid(),
                 startticks=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                 code=code, started_at=datetime.now().astimezone().isoformat())
    lock = out/'owner.lock'
    with lock.open('x', encoding='utf8') as stream:
        json.dump(owner, stream)
    start = time.monotonic()

    def progress(stage, **more):
        status = dict(stage=stage, **owner, time=datetime.now().astimezone().isoformat(),
                      elapsed_seconds=time.monotonic()-start, **more)
        atomic_json(out/'status.json', status)
        print(json.dumps(status), flush=True)

    def deadline(*_):
        raise TimeoutError('Protective worker deadline reached; durable raw blocks preserved')

    signal.signal(signal.SIGALRM, deadline)
    signal.alarm(config['max_seconds'])
    atomic_json(out/'config.resolved.json', config)
    atomic_json(out/'original_method_contract.json', dict(contract=CONTRACT, reference_commit=REFERENCE_COMMIT,
                original_module_sha256=proof, actual_io_adapter_code=code,
                untestable=['original_guide_assignment', 'unfiltered_library_dropout', 'independent_biology'],
                changed_scientific_parameters=False, counts_layer='layers/counts', dense_X_read=False))
    try:
        with threadpool_limits(2), h5py.File(source, 'r') as h:
            progress('read_count_metadata')
            all_obs, var, overlay = adapt_metadata(h, config['catalogue'])
            source_key = dict(path=str(source), bytes=before['bytes'], mtime_ns=before['mtime_ns'],
                              historical_sha256=config['historical_sha256'])
            all_obs.to_csv(out/'tables/source_cell_metadata.csv.gz')
            var.to_csv(out/'tables/all_RNA_identities.csv.gz')
            manifest = cache_counts(h['layers/counts'], out/'raw_count_cache', source_key,
                                    all_obs.index, var.index, progress, config['block_entries'])
        with threadpool_limits(2):
            progress('original_RNA_QC_from_raw_cache')
            mt_mask = var.gene_name.str.upper().str.startswith('MT-').to_numpy()
            totals, detected, mt, _ = moments(out/'raw_count_cache', manifest, len(all_obs), mt_mask)
            if (totals <= 0).any():
                raise ValueError('Empty RNA record')
            metadata_agreement = dict(total_UMI_exact=bool(np.array_equal(totals, all_obs.n_UMI_counts)),
                                      detected_RNA_exact=bool(np.array_equal(detected, all_obs.n_genes)),
                                      UMI_max_absolute_difference=float(np.max(np.abs(totals-all_obs.n_UMI_counts))),
                                      recomputed_from_current_raw_counts=True)
            all_obs['total_counts'] = totals
            all_obs['n_genes_by_counts'] = detected
            all_obs['pct_counts_mt'] = mt/totals*100
            all_obs, bounds = choose_qc(all_obs, {})
            all_obs.to_csv(out/'tables/all_released_cells_with_additional_QC.csv.gz')
            atomic_json(out/'QC_bounds.json', bounds)
            keep = all_obs.main_analysis_stage1.to_numpy(bool)
            obs = all_obs.loc[keep].copy()
            reference = obs.gene_target.eq('NTC').to_numpy()
            if reference.sum() < 2:
                raise ValueError('Too few QC NTC records')
            for field, mask in [('single_guide_before_qc', np.ones(len(all_obs), bool)),
                                ('single_guide_after_stage1_qc', keep)]:
                overlay[field] = overlay.guide.map(all_obs.loc[mask, 'guide_id'].value_counts()).fillna(0).astype(int)
            progress('original_RNA_feature_rule', analysis_cells=len(obs), reference_cells=int(reference.sum()))
            _, _, _, sums = moments(out/'raw_count_cache', manifest, len(all_obs), mt_mask, keep)
            features = select_features(CountMeans(sums, len(obs)), var,
                                       dict(mode='pooled_raw_mean', raw_mean_gt=.5))
            var['raw_mean_after_original_QC'] = sums/len(obs)
            var['response_feature'] = np.isin(np.arange(len(var)), features)
            var.to_csv(out/'tables/all_RNA_feature_selection.csv.gz')
            selected = selected_matrix(out/'raw_count_cache', manifest, keep, features)
            scale = float(obs.loc[reference, 'total_counts'].median())
            progress('original_full_guide_and_10_cell_effects', analysis_cells=len(obs),
                     rna_features=len(features), selected_nnz=int(selected.nnz),
                     normalization_target_umi=scale, guide_count=len(overlay))
            normalized = (sparse.diags(scale/obs.total_counts.to_numpy(dtype=float)) @ selected.astype(float)).tocsr()
            del selected
            names = var.gene_name.to_numpy(str)
            records, means, centers, effects, first, metrics, selections = effect_profiles(
                normalized, obs, overlay, reference, names[features], names,
                min_cells=10, matched_cells=10, repeats=8, seed=20261022, pseudocount=1.)
            records.to_csv(out/'tables/guide_effect_metadata.csv', index=False)
            metrics.to_csv(out/'tables/matched_effect_metrics.csv', index=False)
            atomic_json(out/'matched_cell_selections.json', selections)
            np.savez_compressed(out/'tables/rna_effect_profiles.npz', means=means,
                                reference_means=centers, log2fc=effects, matched_first=first,
                                guides=records.guide.to_numpy(str), gene_ids=var.index.to_numpy(str)[features],
                                feature_indices=features, pseudocount=1.)
            overlay.to_csv(out/'tables/control_annotations.csv', index=False)
            obs.to_csv(out/'tables/analysis_cells.csv.gz')
            var.iloc[features].to_csv(out/'tables/response_features.csv')
            flow = qc_flow(all_obs)
            flow.to_csv(out/'tables/technical_cell_flow.csv', index=False)
            progress('original_target_aggregation', eligible_guides=len(records))
            support, pairs = target_support(normalized, obs, overlay, records, effects, reference)
            support.to_csv(out/'tables/target_support.csv', index=False)
            pairs.to_csv(out/'tables/within_target_pairs.csv', index=False)
            progress('method_discrimination_backgrounds', within_target_pairs=len(pairs))
            backgrounds, coverage = diagnostic_backgrounds(obs, records, effects, pairs, 20261022)
            backgrounds.to_csv(out/'tables/diagnostic_matched_backgrounds.csv', index=False)
            coverage.to_csv(out/'tables/guide_batch_coverage.csv', index=False)
            import matplotlib
            matplotlib.use('Agg')
            figure_annotations()
            progress('original_corresponding_figures')
            original_plots.technical_figures(all_obs, overlay, flow, out)
            examples, example_selections = original_run.example_means(normalized, obs, overlay, reference, 10, 20261022)
            atomic_json(out/'example_cell_selections.json', example_selections)
            original_plots.effect_background_figures(records, effects, first, metrics, examples, out)
            original_run.agreement_figure(support, pairs, out)
            summary = dict(code=code, original_reference_commit=REFERENCE_COMMIT,
                           source_identity=source_key, metadata_agreement=metadata_agreement,
                           source_cells=len(all_obs), additional_RNA_QC_cells=len(obs),
                           response_features=len(features), eligible_guides=len(records),
                           reference_cells=int(reference.sum()), normalization_target_umi=scale,
                           reference_guides=int(obs.loc[reference, 'guide_id'].nunique()),
                           targets_assessable=int(support.assessment.eq('descriptive_assessable').sum()),
                           within_target_pairs=len(pairs), contract=CONTRACT,
                           scope='method_challenge_on_preselected_Strong_release',
                           entire_atlas_quality_assessed=False, functional_hits_called=False,
                           elapsed_seconds=time.monotonic()-start,
                           versions={p: md.version(p) for p in ['numpy','pandas','scipy','h5py','matplotlib']})
            atomic_json(out/'summary.json', summary)
        if identity(source) != before:
            raise ValueError('Source identity changed during method challenge')
        if sha256(config['catalogue']) != config['catalogue_sha256']:
            raise ValueError('Author catalogue changed during method challenge')
        progress('manifest')
        files = {p.relative_to(out).as_posix(): dict(bytes=p.stat().st_size, sha256=sha256(p))
                 for p in sorted(out.rglob('*')) if p.is_file()
                 and p.name not in ['owner.lock','status.json','output_manifest.json']
                 and 'raw_count_cache' not in p.parts}
        atomic_json(out/'output_manifest.json', dict(outputs=files, raw_cache_manifest=sha256(out/'raw_count_cache/manifest.json')))
        progress('complete_pending_method_interpretation', files=len(files))
        return summary
    except BaseException as exc:
        progress('failed', error=repr(exc))
        (out/'traceback.txt').write_text(traceback.format_exc(), encoding='utf8')
        raise
    finally:
        signal.alarm(0)
        lock.unlink()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    run(parser.parse_args().config)
