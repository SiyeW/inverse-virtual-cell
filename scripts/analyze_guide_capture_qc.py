"""Bounded-memory coverage and stage-1 QC for 10x CRISPR Guide Capture.

Read the CSC feature-by-cell 10x H5 in cell blocks. Preserve the full cell
population and raw RNA counts; do not use singleton-only inputs as the QC
denominator. This module does not infer genomic DNA editing efficiency.
"""
from __future__ import annotations

from datetime import datetime
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def file_hash(path: Path, algorithm: str = 'sha256') -> str:
    digest = hashlib.new(algorithm)
    with path.open('rb') as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def feature_table(matrix: h5py.Group) -> pd.DataFrame:
    group = matrix['features']
    return pd.DataFrame({
        'gene_id': group['id'].asstr()[:],
        'gene_name': group['name'].asstr()[:],
        'feature_type': group['feature_type'].asstr()[:],
    })


def iter_cells(matrix: h5py.Group, block_size: int = 128):
    pointers = matrix['indptr'][:]
    cells = len(pointers) - 1
    for start in range(0, cells, block_size):
        stop = min(start + block_size, cells)
        lo, hi = int(pointers[start]), int(pointers[stop])
        data = matrix['data'][lo:hi]
        indices = matrix['indices'][lo:hi]
        if not np.isfinite(data).all() or (data < 0).any() or not np.equal(data, np.floor(data)).all():
            raise ValueError('RNA/guide source contains non-finite, negative, or fractional counts')
        if len(data) and data.max() > np.iinfo(np.int32).max:
            raise OverflowError('Source count exceeds int32')
        if len(indices) and (indices.min() < 0 or indices.max() >= matrix['shape'][0]):
            raise ValueError('Source feature index outside matrix dimensions')
        yield start, stop, data.astype(np.int32, copy=False), indices, pointers[start:stop + 1] - lo


def measure_cells(matrix_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    with h5py.File(matrix_path, 'r') as handle:
        matrix = handle['matrix']
        features = feature_table(matrix)
        barcodes = matrix['barcodes'].asstr()[:]
        if len(set(barcodes)) != len(barcodes):
            raise ValueError('Matrix cell barcodes are duplicated')
        rna = features['feature_type'].eq('Gene Expression').to_numpy()
        guides = features['feature_type'].eq('CRISPR Guide Capture').to_numpy()
        mitochondrial = rna & features['gene_name'].str.upper().str.startswith('MT-').to_numpy()
        if not rna.any() or not mitochondrial.any():
            raise ValueError('RNA or mitochondrial features absent; inspect reference annotation')
        metrics = np.zeros((len(barcodes), 5), dtype=np.int64)
        gene_detected = np.zeros(len(features), dtype=np.int64)
        for start, stop, values, rows, ptr in iter_cells(matrix):
            for local, cell in enumerate(range(start, stop)):
                idx = rows[ptr[local]:ptr[local + 1]]
                val = values[ptr[local]:ptr[local + 1]]
                positive = val > 0
                if len(np.unique(idx)) != len(idx):
                    raise ValueError('Source cell contains duplicated feature indices')
                metrics[cell] = [val[rna[idx]].sum(dtype=np.int64), (positive & rna[idx]).sum(),
                                 val[mitochondrial[idx]].sum(dtype=np.int64), val[guides[idx]].sum(dtype=np.int64),
                                 (positive & rna[idx]).sum()]
                gene_detected[idx[positive & rna[idx]]] += 1
        obs = pd.DataFrame(metrics, index=pd.Index(barcodes, name='cell_barcode'),
                           columns=['total_counts', 'n_genes_by_counts', 'mt_counts', 'raw_guide_umis', 'rna_nnz'])
        obs['pct_counts_mt'] = np.divide(obs['mt_counts'], obs['total_counts'],
                                          out=np.zeros(len(obs), dtype=float), where=obs['total_counts'].to_numpy() > 0) * 100
        features['n_cells_by_counts'] = gene_detected
    return obs, features, rna


def join_guides(obs: pd.DataFrame, calls: pd.DataFrame, reference: pd.DataFrame, control_target: str) -> pd.DataFrame:
    required_calls = {'cell_barcode', 'num_features', 'feature_call', 'num_umis'}
    if not required_calls.issubset(calls.columns) or not {'name', 'target_gene_name'}.issubset(reference.columns):
        raise ValueError('Missing guide-call or reference fields')
    if calls['cell_barcode'].duplicated().any() or reference['name'].duplicated().any():
        raise ValueError('Duplicate barcode or guide reference identity')
    if not set(calls['cell_barcode']).issubset(set(obs.index)):
        raise ValueError('Guide calls contain cells outside the selected expression matrix')
    if reference[['name', 'target_gene_name']].isna().any().any():
        raise ValueError('Missing guide or target names in reference')
    result = obs.copy()
    joined = calls.set_index('cell_barcode').reindex(obs.index)
    result['guide_num_features'] = joined['num_features'].fillna(0).astype(int)
    result['guide_call'] = joined['feature_call'].fillna('').astype(object)
    result['guide_umi_call'] = joined['num_umis'].fillna('').astype(str).astype(object)
    mapping = reference.set_index('name')['target_gene_name'].to_dict()
    targets, singleton_umis = [], []
    for row in result.itertuples():
        names = str(row.guide_call).split('|') if row.guide_call else []
        if len(names) != row.guide_num_features:
            raise ValueError(f'Guide number disagrees with assignment: {row.Index}')
        if any(name not in mapping for name in names):
            raise ValueError(f'Guide absent from authoritative reference: {row.Index}')
        targets.append('|'.join(mapping[name] for name in names))
        singleton_umis.append(int(float(row.guide_umi_call)) if len(names) == 1 else 0)
    result['target_call'] = np.asarray(targets, dtype=object)
    result['guide_status'] = np.select([result['guide_num_features'].eq(1), result['guide_num_features'].gt(1)],
                                      ['single_guide', 'multi_guide'], default='no_confident_call')
    result['guide_id'] = np.where(result['guide_num_features'].eq(1), result['guide_call'], '')
    result['perturbed_gene'] = np.where(result['guide_num_features'].eq(1), result['target_call'], '')
    result.loc[result['perturbed_gene'].eq(control_target), 'perturbed_gene'] = 'ctrl'
    result['guide_umis'] = np.asarray(singleton_umis, dtype=np.int64)
    result['is_control'] = result['perturbed_gene'].eq('ctrl')
    result['guide_detection'] = np.where(result['raw_guide_umis'].eq(0), 'no_guide_molecules', 'guide_signal_detected')
    return result


def mad(values: np.ndarray) -> float:
    return float(np.median(np.abs(values - np.median(values))))


def choose_qc(obs: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, dict]:
    result = obs.copy()
    multiplier = float(config.get('count_mad_multiplier', 5))
    bounds = {}
    for name, floor in [('total_counts', 500), ('n_genes_by_counts', 200)]:
        logged = np.log1p(obs[name].to_numpy())
        center, deviation = float(np.median(logged)), mad(logged)
        bounds[name] = {'low': max(floor, int(np.floor(np.expm1(center - multiplier * deviation) + 1e-8))),
                        'high_diagnostic': float(np.expm1(center + multiplier * deviation)),
                        'log_median': center, 'log_mad': deviation}
        result['low_' + name] = obs[name] < bounds[name]['low']
        result['high_' + name] = obs[name] > bounds[name]['high_diagnostic']
    mt = obs['pct_counts_mt'].to_numpy()
    mt_limit = max(float(config.get('mitochondrial_floor_pct', 8)), float(np.median(mt) + 3 * mad(mt)))
    bounds['pct_counts_mt'] = {'high': mt_limit, 'median': float(np.median(mt)), 'mad': mad(mt)}
    result['high_mitochondrial'] = obs['pct_counts_mt'] > mt_limit
    result['qc_pass_stage1'] = ~(result['low_total_counts'] | result['low_n_genes_by_counts'] | result['high_mitochondrial'])
    result['main_analysis_stage1'] = result['qc_pass_stage1'] & result['guide_status'].eq('single_guide')
    reasons = []
    for row in result.itertuples():
        reasons.append(';'.join(label for flag, label in [(row.low_total_counts, 'low_counts'),
                         (row.low_n_genes_by_counts, 'low_genes'), (row.high_mitochondrial, 'high_mt')] if flag))
    result['qc_exclusion_reason'] = np.asarray(reasons, dtype=object)
    result['doublet_assessment'] = 'pending'
    bounds['policy'] = 'Low RNA counts/genes below lower 5 MAD(log1p), with floors 500/200; mt above both 3 MAD and 8%; high RNA diagnostics retained pending doublet assessment.'
    return result, bounds


def coverage_tables(obs: pd.DataFrame, reference: pd.DataFrame, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    output.mkdir(parents=True, exist_ok=True)
    guide = reference.copy().set_index('name')
    single = obs['guide_status'].eq('single_guide')
    masks = {'single_guide_before_qc': single, 'single_guide_after_stage1_qc': obs['main_analysis_stage1']}
    for label, mask in masks.items():
        guide[label] = obs.loc[mask, 'guide_id'].value_counts().reindex(guide.index, fill_value=0).astype(int)
    assigned = [name for text in obs.loc[obs['guide_num_features'].gt(0), 'guide_call'] for name in text.split('|')]
    guide['cells_with_guide_including_multi'] = pd.Series(assigned, dtype=object).value_counts().reindex(guide.index, fill_value=0).astype(int)
    guide['guides_target'] = guide['target_gene_name']
    targets = pd.DataFrame(index=pd.Index(reference['target_gene_name'].drop_duplicates(), name='target_gene_name'))
    targets['library_guides'] = reference.groupby('target_gene_name').size().reindex(targets.index)
    for label in masks:
        targets[label] = guide.groupby('target_gene_name')[label].sum().reindex(targets.index).astype(int)
    targets['observed_guides_after_stage1_qc'] = guide[guide['single_guide_after_stage1_qc'].gt(0)].groupby('target_gene_name').size().reindex(targets.index, fill_value=0)
    targets['is_control'] = targets.index == 'Non-Targeting'
    targets['at_least_20_cells'] = targets['single_guide_after_stage1_qc'].ge(20)
    targets['at_least_50_cells'] = targets['single_guide_after_stage1_qc'].ge(50)
    targets['at_least_100_cells'] = targets['single_guide_after_stage1_qc'].ge(100)
    guide.to_csv(output / 'coverage_per_sgrna.csv', encoding='utf-8-sig')
    targets.to_csv(output / 'coverage_per_target.csv', encoding='utf-8-sig')
    pd.crosstab(obs['guide_status'], obs['qc_pass_stage1']).to_csv(output / 'guide_assignment_by_qc.csv', encoding='utf-8-sig')
    obs.to_csv(output / 'cell_qc_and_guide_assignments.csv.gz', compression='gzip')
    return guide, targets


def write_streamed_h5ad(matrix_path: Path, obs: pd.DataFrame, features: pd.DataFrame,
                       rna: np.ndarray, selected: np.ndarray, path: Path, provenance: dict) -> dict:
    """Write valid AnnData CSR counts in bounded memory, including a counts layer."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.partial.h5ad')
    selected = np.asarray(selected, dtype=bool)
    if len(selected) != len(obs) or not selected.any():
        raise ValueError('Empty or invalid cell selection')
    variables = features.loc[rna].copy()
    variables.index = pd.Index(variables['gene_id'].astype(str), name='feature_id')
    # Gene IDs, including custom PB identifiers, remain the primary feature key.
    if not variables.index.is_unique:
        raise ValueError('RNA feature IDs are duplicated')
    frame = obs.loc[selected].copy()
    for name in frame.select_dtypes(include=['object', 'string']).columns:
        frame[name] = pd.Categorical(frame[name])
    for name in variables.select_dtypes(include=['object', 'string']).columns:
        variables[name] = pd.Categorical(variables[name])
    n_cells, n_genes = len(frame), len(variables)
    n_entries = int(obs.loc[selected, 'rna_nnz'].sum())
    template = ad.AnnData(X=sparse.csr_matrix((n_cells, n_genes), dtype=np.int32), obs=frame, var=variables)
    template.uns['analysis_provenance'] = provenance
    template.write_h5ad(temporary)
    del template
    mapping = np.full(len(rna), -1, dtype=np.int32)
    mapping[rna] = np.arange(n_genes, dtype=np.int32)
    offset, row_number = 0, 0
    ptr = np.zeros(n_cells + 1, dtype=np.int64)
    with h5py.File(matrix_path, 'r') as source, h5py.File(temporary, 'r+') as target:
        x = target['X']
        for name in ['data', 'indices', 'indptr']:
            del x[name]
        options = {'shape': (n_entries,), 'chunks': (min(262144, max(n_entries, 1)),), 'compression': 'lzf'}
        values_out = x.create_dataset('data', dtype=np.int32, **options)
        indices_out = x.create_dataset('indices', dtype=np.int32, **options)
        for start, stop, values, indices, pointers in iter_cells(source['matrix']):
            buffers_values, buffers_indices = [], []
            for local, cell in enumerate(range(start, stop)):
                if not selected[cell]:
                    continue
                idx = indices[pointers[local]:pointers[local + 1]]
                val = values[pointers[local]:pointers[local + 1]]
                keep = rna[idx] & (val > 0)
                mapped = mapping[idx[keep]]
                counts = val[keep]
                order = np.argsort(mapped)
                buffers_values.append(counts[order])
                buffers_indices.append(mapped[order])
                ptr[row_number + 1] = ptr[row_number] + len(mapped)
                row_number += 1
            if buffers_values:
                data = np.concatenate(buffers_values)
                rows = np.concatenate(buffers_indices)
                values_out[offset:offset + len(data)] = data
                indices_out[offset:offset + len(rows)] = rows
                offset += len(data)
        if offset != n_entries or row_number != n_cells:
            raise ValueError('Written CSR dimensions do not match selected cell metrics')
        x.create_dataset('indptr', data=ptr)
        # Same immutable count object; HDF5 hard link avoids duplicating disk bytes.
        target['layers']['counts'] = x
        target.flush()
    with h5py.File(temporary, 'r') as verify:
        x = verify['X']
        if tuple(x.attrs['shape']) != (n_cells, n_genes) or x['indptr'][-1] != n_entries:
            raise ValueError('Reopened H5AD dimensions are inconsistent')
        if verify['layers/counts'].id != x.id or x['data'].dtype != np.dtype('int32'):
            raise ValueError('Count-layer encoding failed')
        for start in range(0, n_cells, 128):
            stop = min(start + 128, n_cells)
            lo, hi = int(ptr[start]), int(ptr[stop])
            block = sparse.csr_matrix((x['data'][lo:hi], x['indices'][lo:hi], ptr[start:stop + 1] - lo), shape=(stop - start, n_genes))
            expected = frame['total_counts'].iloc[start:stop].to_numpy()
            if not np.array_equal(np.asarray(block.sum(axis=1)).ravel(), expected):
                raise ValueError('Reopened H5AD row counts differ from source QC metrics')
    temporary.replace(path)
    return {'path': str(path), 'cells': n_cells, 'features': n_genes, 'nnz': n_entries,
            'bytes': path.stat().st_size, 'sha256': file_hash(path), 'reopened_row_sums_match': True}


def plot_qc(obs: pd.DataFrame, guide: pd.DataFrame, targets: pd.DataFrame, thresholds: dict, output: Path) -> list[str]:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    colors = {'single_guide': '#367DA2', 'multi_guide': '#C87840', 'no_confident_call': '#8B9399'}
    for status, cells in obs.groupby('guide_status', observed=True):
        axes[0, 0].scatter(cells['total_counts'], cells['n_genes_by_counts'], s=2, alpha=.25, c=colors[status], label=status, rasterized=True)
    axes[0, 0].set(xscale='log', yscale='log', xlabel='RNA UMIs', ylabel='Detected RNA features')
    axes[0, 0].legend(markerscale=3, fontsize=8)
    for axis, metric in [(axes[0, 1], 'total_counts'), (axes[0, 2], 'n_genes_by_counts')]:
        axis.hist(np.log10(1 + obs[metric]), bins=70, color='#367DA2')
        axis.axvline(np.log10(1 + thresholds[metric]['low']), color='#A33A3A', linestyle='--')
        axis.set(xlabel='log10(1 + ' + metric + ')', ylabel='Cells')
    axes[1, 0].hist(obs['pct_counts_mt'], bins=70, color='#367DA2')
    axes[1, 0].axvline(thresholds['pct_counts_mt']['high'], color='#A33A3A', linestyle='--')
    axes[1, 0].set(xlabel='Mitochondrial RNA UMIs (%)', ylabel='Cells')
    counts = obs['guide_status'].value_counts().reindex(list(colors), fill_value=0)
    axes[1, 1].bar(range(len(counts)), counts, color=list(colors.values()))
    axes[1, 1].set_xticks(range(len(counts)), ['Single guide', 'Multi guide', 'No confident call'], rotation=15)
    for index, number in enumerate(counts):
        axes[1, 1].text(index, number, f'{number:,}', ha='center', va='bottom')
    stages = [len(obs), int(obs['qc_pass_stage1'].sum()), int(obs['main_analysis_stage1'].sum())]
    axes[1, 2].bar(['All', 'RNA QC pass', 'Single + QC'], stages, color=['#8B9399', '#367DA2', '#418C70'])
    axes[1, 2].set(ylabel='Cells')
    fig.suptitle('Chip 01: coverage and stage-1 QC (doublet assessment pending)')
    path = output / '01_cell_qc_overview.png'
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(str(path))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for axis, frame, title in [(axes[0], guide, 'Per sgRNA (all library guides)'),
                                (axes[1], targets.loc[~targets['is_control']], 'Per target (excluding NTC)')]:
        a = np.sort(frame['single_guide_before_qc'].to_numpy())[::-1]
        b = np.sort(frame['single_guide_after_stage1_qc'].to_numpy())[::-1]
        axis.plot(np.arange(1, len(a) + 1), a, label='Before RNA QC', c='#8B9399')
        axis.plot(np.arange(1, len(b) + 1), b, label='After stage-1 RNA QC', c='#367DA2')
        axis.set(xlabel='Coverage rank (independently sorted)', ylabel='Single-guide cells', yscale='symlog', title=title)
        axis.legend()
    path = output / '02_perturbation_coverage.png'
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(str(path))
    return paths


def run_analysis(config_path: Path) -> dict:
    import psutil
    config = json.loads(config_path.read_text(encoding='utf-8'))
    root = Path(config['output_root'])
    for name in ['data', 'tables', 'figures', 'reports', 'provenance']:
        (root / name).mkdir(parents=True, exist_ok=True)
    status_path = root / 'provenance/analysis_status.json'
    def stage(name):
        write_json(status_path, {'stage': name, 'time': datetime.now().astimezone().isoformat(),
                                 'rss_bytes': psutil.Process().memory_info().rss})
        print(f'{datetime.now().isoformat(timespec="seconds")} {name}', flush=True)
    stage('measure_full_matrix')
    matrix_path = Path(config['matrix_h5'])
    obs, features, rna = measure_cells(matrix_path)
    calls = pd.read_csv(config['protospacer_calls'])
    reference = pd.read_csv(config['feature_reference'])
    obs = join_guides(obs, calls, reference, config.get('control_target', 'Non-Targeting'))
    obs['chip_id'] = config.get('chip_id', 'chip01')
    stage('choose_and_record_stage1_qc')
    obs, thresholds = choose_qc(obs, config)
    write_json(root / 'provenance/qc_thresholds.json', thresholds)
    guide, targets = coverage_tables(obs, reference, root / 'tables')
    pd.DataFrame([{'metric': metric, 'median': float(obs[metric].median()),
                   'p01': float(obs[metric].quantile(.01)), 'p99': float(obs[metric].quantile(.99))}
                  for metric in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'raw_guide_umis']]).to_csv(root / 'tables/cell_metric_summary.csv', index=False)
    company = pd.read_csv(config['metrics_summary'])
    company.to_csv(root / 'tables/company_metrics_summary.csv', index=False, encoding='utf-8-sig')
    feature_types = features.groupby('feature_type').size().astype(int).to_dict()
    versions = {p: importlib.metadata.version(p) for p in ['numpy', 'scipy', 'pandas', 'anndata', 'h5py']}
    code_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=Path(__file__).resolve().parents[1], text=True).strip()
    provenance = {'chip_id': config.get('chip_id', 'chip01'), 'code_commit': code_commit,
                  'stage': 'coverage_and_stage1_qc', 'doublet_assessment': 'pending',
                  'DNA_editing_efficiency': 'not_assessed_from_RNA_guide_assay',
                  'counts_encoding': 'raw int32 CSR; counts layer is an HDF5 hard link to X',
                  'input_hashes_json': str(config['input_hashes_json'])}
    outputs = []
    for name, selected in [('chip01_all_cells_stage1_qc.h5ad', np.ones(len(obs), dtype=bool)),
                            ('chip01_single_guide_stage1_qc.h5ad', obs['main_analysis_stage1'].to_numpy())]:
        stage('write_' + name)
        outputs.append(write_streamed_h5ad(matrix_path, obs, features, rna, selected, root / 'data' / name, provenance))
    stage('plot_qc_and_coverage')
    figures = plot_qc(obs, guide, targets, thresholds, root / 'figures')
    noncontrols = targets.loc[~targets['is_control']]
    summary = {'cells': len(obs), 'rna_features': int(rna.sum()), 'feature_types': feature_types,
               'guide_status': {str(k): int(v) for k, v in obs['guide_status'].value_counts().items()},
               'rna_qc_pass_cells': int(obs['qc_pass_stage1'].sum()),
               'main_stage1_cells': int(obs['main_analysis_stage1'].sum()),
               'controls_before_qc': int(obs['is_control'].sum()),
               'controls_after_qc': int((obs['is_control'] & obs['main_analysis_stage1']).sum()),
               'library_guides': len(reference), 'library_targets_including_control': len(targets),
               'observed_guides_after_qc': int(guide['single_guide_after_stage1_qc'].gt(0).sum()),
               'noncontrol_targets_with_20_cells': int(noncontrols['at_least_20_cells'].sum()),
               'noncontrol_targets_with_50_cells': int(noncontrols['at_least_50_cells'].sum()),
               'noncontrol_targets_with_100_cells': int(noncontrols['at_least_100_cells'].sum()),
               'no_guide_molecules_from_matrix': int(obs['raw_guide_umis'].eq(0).sum()),
               'doublet_assessment': 'pending', 'outputs': outputs, 'figures': figures,
               'thresholds': thresholds, 'versions': versions, 'code_commit': code_commit,
               'rss_bytes_at_finish': psutil.Process().memory_info().rss}
    write_json(root / 'provenance/summary.json', summary)
    report = f'''# 第 1 张芯片：覆盖度与第一阶段质控

全量 {len(obs):,} 个细胞，{int(rna.sum()):,} 个 RNA feature（含公司自定义参考，不全部等同于标准基因）。
单 guide {summary['guide_status'].get('single_guide', 0):,}；多 guide {summary['guide_status'].get('multi_guide', 0):,}；无可靠分配 {summary['guide_status'].get('no_confident_call', 0):,}。

## 质控与覆盖度

第一阶段 RNA 质控通过 {summary['rna_qc_pass_cells']:,} 个细胞，单 guide 主分析候选 {summary['main_stage1_cells']:,} 个，其中阴性对照 {summary['controls_after_qc']:,} 个。
文库有 {len(reference):,} 条 guide、{len(targets):,} 个目标标签（含阴性对照）。质控后有至少 20/50/100 个单 guide 细胞的非对照目标分别为 {summary['noncontrol_targets_with_20_cells']}/{summary['noncontrol_targets_with_50_cells']}/{summary['noncontrol_targets_with_100_cells']}。
20/50/100 只是覆盖度分层，不能当作统计功效保证。

低计数阈值 {thresholds['total_counts']['low']:.1f} UMI；低 feature 阈值 {thresholds['n_genes_by_counts']['low']:.1f}；线粒体比例上限 {thresholds['pct_counts_mt']['high']:.2f}%。完整阈值依据、逐细胞排除原因及过滤前后覆盖度见 provenance 和 tables。

## 结果边界与下一步

本轮完成全量 guide 分配、覆盖度、常规 RNA 质控、稀疏原始计数 h5ad 及图表。高计数/高 feature 细胞只标记诊断，不自动当成双细胞排除。
双细胞评分和表达降维尚未完成；主分析候选为 stage-1 版本，不宣称已经完成最终双细胞质控。
后续依据实际资源完成双细胞诊断、表达结构检查，再开展能量距离、差异表达和功能分析。
单张芯片不自动提供独立生物学重复。DNA indel/frameshift 不能从本 RNA/guide 数据可靠评估，记为未评估。

## 可追溯性

代码提交：`{code_commit}`。输入哈希见 `{config['input_hashes_json']}`。输出逐行计数已重新读取核对，SHA-256 见 summary.json。
QC 策略参考：https://www.sc-best-practices.org/preprocessing-visualization/quality-control/ 。
'''
    (root / 'reports/01_coverage_and_stage1_qc.md').write_text(report, encoding='utf-8')
    stage('complete_stage1')
    return summary
