"""One CPU-only entry point: adopted descriptive Perturb-seq QC, figures and report.

Usage: run_core_analysis.py run --config local.json [--smoke]
       run_core_analysis.py verify --output /path/to/completed/run
       run_core_analysis.py prepare --config local-10x.json
"""
import argparse
from datetime import datetime
import importlib.metadata as md
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

from perturbseq_core_io import (atomic_json, identity, load_config, read_catalogue,
    read_source, revision, select_features, sha256)
from perturbseq_core_effects import TYPES, effect_profiles, stable_rng, target_support, qc_flow
from perturbseq_core_plots import (COLORS, NAMES, technical_figures,
    effect_background_figures, target_atlas, save)
from perturbseq_core_report import write_report, verify_output


def role_summary(overlay, records):
    rows = []
    for role in TYPES:
        group = overlay.loc[overlay.source_control_type.eq(role)]
        values = group.single_guide_after_stage1_qc.to_numpy()
        rows.append(dict(source_control_type=role, library_guides=len(group), cells=int(values.sum()),
            eligible_guides=int(records.source_control_type.eq(role).sum()),
            zero_guides=int((values == 0).sum()), median_cells=float(np.median(values)) if len(values) else None))
    return pd.DataFrame(rows)


def example_means(norm, obs, overlay, reference, matched, seed):
    guides = obs.guide_id.astype(str).to_numpy()
    examples, selections = [], []
    groups = []
    for role in TYPES[:2]:
        members = overlay.loc[overlay.source_control_type.eq(role)].sort_values(
            ['single_guide_after_stage1_qc', 'guide'], ascending=[False, True])
        eligible = members.loc[members.single_guide_after_stage1_qc.ge(matched)]
        if len(eligible):
            row = eligible.iloc[0]
            groups.append((row.guide, role, guides == row.guide, row.guide))
    for target, members in overlay.loc[overlay.source_control_type.eq('positive_control')].groupby('target_gene_name', sort=True):
        mask = np.isin(guides, members.guide)
        if mask.sum() >= matched:
            groups.append((str(target), 'positive_control', mask, None))
    for label, role, mask, own_guide in groups:
        ref = reference & ~mask
        if ref.sum() < 2:
            raise ValueError('Insufficient nonoverlapping example reference')
        # Same deterministic key convention as the accepted report.
        chosen = np.sort(stable_rng(seed, 'example:' + label).choice(np.flatnonzero(mask), matched, replace=False))
        examples.append(dict(label=f'{NAMES[role]}: {label}; n={matched}',
            mean=np.asarray(norm[chosen].mean(0)).ravel(), reference_mean=np.asarray(norm[ref].mean(0)).ravel(), color=COLORS[role]))
        selections.append(dict(group=label, role=role, cells=matched,
            source_rows=obs.iloc[chosen].source_row.astype(int).tolist(), reference_cells=int(ref.sum()), reference_overlap_cells=0))
    return examples, selections


def agreement_figure(support, pairs, output):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
    bins = [int(support.eligible_guides.eq(i).sum()) for i in range(3)] + [int(support.eligible_guides.ge(3).sum())]
    axes[0].bar(['0', '1', '2', '3+'], bins, color='#3c7c9f')
    axes[0].set(xlabel='Eligible sgRNAs per target', ylabel='Targets', title='Target coverage')
    shown = 0
    for role in ['positive_control', 'targeting']:
        part = support.loc[support.source_control_type.eq(role) & support.assessment.eq('descriptive_assessable')]
        if len(part):
            axes[1].scatter(part.largest_guide_fraction, part.median_pair_r, s=24, alpha=.7,
                            color=COLORS[role], label=NAMES[role])
            shown += len(part)
    axes[1].axhline(0, color='grey', lw=.8)
    if shown:
        axes[1].legend(frameon=False, fontsize=8)
    else:
        axes[1].text(.5, .5, 'No assessable targets', ha='center', transform=axes[1].transAxes)
    axes[1].set(xlabel='Largest sgRNA cell fraction', ylabel='Median Pearson r', title='RNA-response agreement')
    save(fig, output/'figures/07_target_support.png')


def prepare(config_path):
    """Wrap the established 10x stage with fixed inputs and exclusive new output."""
    from analyze_guide_capture_qc import run_analysis
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text(encoding='utf8'))
    import re
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}', config['chip_id']):
        raise ValueError('Unsafe chip identifier')
    code = revision(config['expected_commit'])
    root = Path(config['output_root']).resolve()
    repo = Path(__file__).resolve().parents[1]
    if root.parent != Path(config['allowed_output_parent']).resolve() or root.is_relative_to(repo):
        raise ValueError('Distinct external approved output required')
    keys = ['matrix_h5','protospacer_calls','feature_reference','metrics_summary']
    paths = [Path(config[key]).resolve() for key in keys]
    for key,path in zip(keys,paths):
        if not Path(config[key]).is_absolute():
            raise ValueError('Preparation input paths must be absolute: '+key)
        config[key]=str(path)
    config['output_root']=str(root)
    if any(p == root or p.is_relative_to(root) for p in paths):
        raise ValueError('Output must not contain an input')
    if root.exists():
        raise FileExistsError('Existing baseline is protected; inspect it instead of rerunning')
    snapshots = {key:identity(path) for key,path in zip(keys,paths)}
    for key,path in zip(keys,paths):
        if sha256(path) != config['input_sha256'][key]:
            raise ValueError('10x source SHA256 differs: '+key)
    root.mkdir()
    lock=root/'owner.lock'
    with lock.open('x',encoding='utf8') as stream:
        json.dump(dict(pid=os.getpid(),code=code),stream)
    try:
        atomic_json(root/'config.resolved.json',config)
        atomic_json(root/'status.json',dict(stage='preparing',phase='prepare',code=code))
        with threadpool_limits(config.get('threads',2)):
            summary=run_analysis(root/'config.resolved.json')
        for key,path in zip(keys,paths):
            if identity(path)!=snapshots[key]:
                raise ValueError('10x source identity changed: '+key)
        atomic_json(root/'summary.json',dict(summary,phase='prepare',code=code,source_modified=False))
        # Large H5ADs were already reopened, row-count-checked and hashed by the writer.
        proofs={Path(p['path']).resolve():p for p in summary['outputs']}
        manifest={}
        for path in sorted(root.rglob('*')):
            if not path.is_file() or path.name in ['owner.lock','status.json','output_manifest.json']:continue
            prior=proofs.get(path.resolve())
            manifest[path.relative_to(root).as_posix()]=dict(bytes=path.stat().st_size,
                sha256=prior['sha256'] if prior else sha256(path))
        atomic_json(root/'output_manifest.json',manifest)
        atomic_json(root/'status.json',dict(stage='complete',phase='prepare',code=code))
    except BaseException as error:
        atomic_json(root/'status.json',dict(stage='failed',phase='prepare',error=repr(error)))
        raise
    finally:
        lock.unlink()
    return summary


def run(config_path, smoke=False):
    config = load_config(config_path)
    code = revision(config['expected_commit'], smoke)
    out = Path(config['output'])
    repo = Path(__file__).resolve().parents[1]
    if out.parent != Path(config['allowed_output_parent']) or out.is_relative_to(repo):
        raise ValueError('A distinct external output directly under allowed_output_parent is required')
    sources = [Path(config[key]['path']) for key in ['source', 'catalogue']]
    if config['features']['mode'] == 'frozen_ids':
        sources.append(Path(config['features']['path']))
    if any(p == out or p.is_relative_to(out) for p in sources):
        raise ValueError('Output must not contain or replace any input')
    for source in sources:
        if not source.is_file():
            raise ValueError('Missing input: ' + str(source))
    out.mkdir(exist_ok=False)  # An existing/failed run is never silently resumed or replaced.
    lock = out/'owner.lock'
    owner = dict(pid=os.getpid(), started_at=datetime.now().astimezone().isoformat(), code=code)
    with lock.open('x', encoding='utf8') as stream:
        json.dump(owner, stream)
    started = time.monotonic()

    def state(stage, **more):
        status = dict(stage=stage, elapsed_seconds=time.monotonic()-started, time=datetime.now().astimezone().isoformat(), **owner, **more)
        atomic_json(out/'status.json', status)
        print(json.dumps(status, ensure_ascii=False), flush=True)

    try:
        import matplotlib
        matplotlib.use('Agg')
        for folder in ['tables', 'figures/targets']:
            (out/folder).mkdir(parents=True)
        atomic_json(out/'config.resolved.json', config)
        state('read_immutable_source')
        with threadpool_limits(config['threads']):
            counts, obs, all_obs, var, source_identity = read_source(config['source'], config['smoke_prefix_rows'] if smoke else None)
            overlay = read_catalogue(config['catalogue'], all_obs, config['control_target'])
            reference = obs.guide_id.astype(str).map(overlay.set_index('guide').source_control_type).eq(TYPES[0]).to_numpy()
            if reference.sum() < 2 or obs.loc[reference, 'guide_id'].nunique() < 2:
                raise ValueError('At least two distinct observed NTC guides needed for whole-guide leave-out')
            scale = config['normalization_target_umi']
            if scale == 'median_ntc':
                scale = float(obs.loc[reference, 'total_counts'].median())
            scale = float(scale)
            if not np.isfinite(scale) or scale <= 0:
                raise ValueError('Positive fixed normalization depth required')
            features = select_features(counts, var, config['features'])
            normalized = (sparse.diags(scale/obs.total_counts.to_numpy(dtype=float)) @ counts[:, features].astype(float)).tocsr()
            names = var.gene_name.astype('string').fillna('').to_numpy(dtype=str)
            state('descriptive_rna_effects', cells=len(obs), rna_features=len(features))
            records, means, centers, effects, first, metrics, selections = effect_profiles(
                normalized, obs, overlay, reference, names[features], names, config['min_guide_cells'],
                config['matched_cells'], config['matched_repeats'], config['seed'])
            overlay.to_csv(out/'tables/control_annotations.csv', index=False)
            obs.to_csv(out/'tables/analysis_cells.csv.gz')
            feature_table=var.iloc[features].copy()
            if 'gene_id' in feature_table:
                if not np.array_equal(feature_table.gene_id.astype(str),feature_table.index.astype(str)):
                    raise ValueError('RNA gene_id column/index identity conflict')
                feature_table=feature_table.drop(columns='gene_id')
            feature_table.to_csv(out/'tables/response_features.csv', index_label='gene_id')
            flow = qc_flow(all_obs)
            flow.to_csv(out/'tables/technical_cell_flow.csv', index=False)
            records.to_csv(out/'tables/guide_effect_metadata.csv', index=False)
            metrics.to_csv(out/'tables/matched_effect_metrics.csv', index=False)
            atomic_json(out/'matched_cell_selections.json', selections)
            np.savez_compressed(out/'tables/rna_effect_profiles.npz', means=means, reference_means=centers,
                log2fc=effects, matched_first=first, guides=records.guide.to_numpy(dtype=str),
                gene_ids=var.index.to_numpy(dtype=str)[features], feature_indices=features, pseudocount=1.)
            roles = role_summary(overlay, records)
            roles.to_csv(out/'tables/control_type_summary.csv', index=False)
            examples, example_selections = example_means(normalized, obs, overlay, reference, config['matched_cells'], config['seed'])
            atomic_json(out/'example_cell_selections.json', example_selections)
            state('figures_and_target_cards')
            technical_figures(all_obs, overlay, flow, out)
            background = effect_background_figures(records, effects, first, metrics, examples, out)
            support, pairs = target_support(normalized, obs, overlay, records, effects, reference)
            support.to_csv(out/'tables/target_support.csv', index=False)
            pairs.to_csv(out/'tables/within_target_pairs.csv', index=False)
            agreement_figure(support, pairs, out)
            cards = target_atlas(overlay, records, effects, support, pairs, out,
                                 max_cards=config['max_target_cards'], min_cells=config['min_guide_cells'])
        if identity(config['source']['path']) != {k: source_identity[k] for k in ['bytes','mtime_ns','device','inode']}:
            raise ValueError('Source identity changed; preserve the run, no automatic retry')
        for spec in [config['catalogue']] + ([config['features']] if config['features']['mode']=='frozen_ids' else []):
            if sha256(spec['path']) != spec['sha256']:
                raise ValueError('Small frozen input changed during analysis')
        summary = dict(run_id=config['run_id'], chip_id=config['chip_id'], code=code,
            scientific_result=not smoke, source_identity=source_identity, counts_modified=False,
            source_cells=len(all_obs), rna_qc_cells=int(all_obs.qc_pass_stage1.sum()), analysis_cells=len(obs),
            normalization_target_umi=scale, source_rna_features=len(var), response_features=len(features),
            feature_policy=config['features'], reference_cells=int(reference.sum()), reference_guides=int(obs.loc[reference,'guide_id'].nunique()),
            eligible_guides=len(records), targets=len(support), target_cards=len(cards),
            targets_assessable=int(support.assessment.eq('descriptive_assessable').sum()) if len(support) else 0,
            matched_cells=config['matched_cells'], matched_repeats=config['matched_repeats'],
            input_scope=config['input_scope'], report_notice=config.get('report_notice', ''),
            max_target_cards=config['max_target_cards'],
            validated_functional_hits=False, inferential_de='not_run', dna_editing_efficiency='not_assessed',
            candidate_positive_control_expected_phenotype='not established by source annotation',
            versions={p:md.version(p) for p in ['numpy','pandas','scipy','anndata','h5py','matplotlib']},
            python=sys.version, launch_command=sys.argv, elapsed_seconds=time.monotonic()-started)
        atomic_json(out/'summary.json', summary)
        write_report(out, summary, roles, support, pairs)
        manifest = {p.relative_to(out).as_posix():dict(bytes=p.stat().st_size, sha256=sha256(p))
                    for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json']}
        atomic_json(out/'output_manifest.json', manifest)
        state('complete', files=len(manifest), scientific_result=not smoke)
    except BaseException as error:
        state('failed', error=repr(error))
        (out/'traceback.txt').write_text(traceback.format_exc(), encoding='utf8')
        raise
    finally:
        lock.unlink()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    execute = commands.add_parser('run', help='New, non-overwriting run from immutable all-cell CSR counts')
    execute.add_argument('--config', type=Path, required=True)
    execute.add_argument('--smoke', action='store_true', help='Contiguous bounded prefix; never a scientific result')
    upstream = commands.add_parser('prepare', help='Verified 10x inputs to all-cell sparse raw-count QC H5AD; not Cell Ranger')
    upstream.add_argument('--config', type=Path, required=True)
    verify = commands.add_parser('verify', help='Read result files only; never rerun analysis or reread raw data')
    verify.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.config, args.smoke)
    elif args.command == 'prepare':
        prepare(args.config)
    else:
        print(json.dumps(verify_output(args.output), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
