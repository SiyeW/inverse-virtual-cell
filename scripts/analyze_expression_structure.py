"""Read-only, bounded-memory expression diagnostics; not final doublet QC.

Uses streaming log-CP10k moments for a documented variable-feature heuristic
and centered sparse PCA. Preserves all source counts and cell identities.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, svds
from threadpoolctl import threadpool_limits


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        while block := stream.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf8')
    tmp.replace(path)


def blocks(handle, rows, n_vars, block_size=128):
    """Yield CSR row blocks without materializing the full source matrix."""
    g = handle['X']
    if g.attrs.get('encoding-type') != 'csr_matrix':
        raise ValueError('Expected CSR counts input')
    for start in range(0, len(rows), block_size):
        picked = rows[start:start + block_size]
        # Read a contiguous source span, then select the requested rows.
        lo_row, hi_row = int(picked[0]), int(picked[-1]) + 1
        ptr = g['indptr'][lo_row:hi_row + 1]
        lo, hi = int(ptr[0]), int(ptr[-1])
        x = sparse.csr_matrix((g['data'][lo:hi], g['indices'][lo:hi], ptr - lo),
                              shape=(hi_row - lo_row, n_vars))
        yield picked, x[picked - lo_row]


def normalized_log(x, totals):
    if np.any(np.asarray(totals) <= 0):
        raise ValueError('Nonpositive RNA total')
    x = x.astype(np.float64)
    x = sparse.diags(10000.0 / np.asarray(totals)) @ x
    x.data = np.log1p(x.data)
    return x.tocsr()


def centered_pca(x, components=30, seed=0):
    """SVD of a centered sparse matrix without a dense cells x genes array."""
    n, p = x.shape
    k = min(components, n - 1, p - 1)
    if k < 1:
        raise ValueError('Too few cells or variable features')
    mean = np.asarray(x.mean(axis=0)).ravel()
    def mv(v):
        return x @ v - np.dot(mean, v)
    def rmv(v):
        return x.T @ v - mean * np.sum(v)
    def mm(v):
        return x @ v - (mean @ v)[None, :]
    def rmm(v):
        return x.T @ v - mean[:, None] * np.sum(v, axis=0)[None, :]
    op = LinearOperator(x.shape, matvec=mv, rmatvec=rmv, matmat=mm,
                        rmatmat=rmm, dtype=np.float64)
    u, singular, vt = svds(op, k=k, random_state=seed)
    order = np.argsort(singular)[::-1]
    singular, u, vt = singular[order], u[:, order], vt[order]
    return u * singular, vt, singular ** 2 / (n - 1)


def run(source, output, hvg_count=1500):
    if hvg_count < 3:
        raise ValueError('hvg_count must be at least 3')
    output.mkdir(parents=True, exist_ok=True)
    status = output / 'expression_diagnostics_status.json'
    if status.exists():
        raise FileExistsError('Existing run status: inspect before choosing a new run directory')
    lock = output / 'expression_diagnostics.lock'
    # Exclusive creation: even if the lock is stale, investigate instead of overriding.
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    def state(stage, **values):
        atomic_json(status, {'stage': stage, 'pid': os.getpid(),
                    'time': datetime.now().astimezone().isoformat(), **values})
    try:
        state('input_verification')
        baseline = source.parent.parent
        summary = json.loads((baseline / 'provenance' / 'summary.json').read_text(encoding='utf8'))
        identity = sha256(source)
        accepted = json.loads((baseline / 'provenance' / 'acceptance_check.json').read_text(encoding='utf8'))
        if not all(x['sha256_matches'] for x in accepted['outputs']):
            raise ValueError('Baseline acceptance failed')
        expected = [entry['sha256'] for entry in summary['outputs']
                    if Path(entry['path']).resolve() == source.resolve()]
        if expected != [identity]:
            raise ValueError('Source path/hash does not match accepted baseline summary')
        a = ad.read_h5ad(source, backed='r')
        obs, var = a.obs.copy(), a.var.copy()
        a.file.close()
        rows = np.flatnonzero(obs['qc_pass_stage1'].to_numpy())
        total = obs['total_counts'].to_numpy()
        standard = var['gene_id'].astype(str).str.match(r'^ENSG\d+(\.\d+)?$').to_numpy()
        excluded_mt = var['gene_name'].astype(str).str.upper().str.startswith('MT-').to_numpy()
        sums, squares, detected = [np.zeros(len(var), dtype=np.float64) for _ in range(3)]
        state('streaming_feature_moments', cells=len(rows))
        with h5py.File(source, 'r') as handle:
            for picked, x in blocks(handle, rows, len(var)):
                z = normalized_log(x, total[picked])
                sums += np.asarray(z.sum(axis=0)).ravel()
                squares += np.asarray(z.multiply(z).sum(axis=0)).ravel()
                detected += np.asarray((x > 0).sum(axis=0)).ravel()
        means = sums / len(rows)
        variance = np.maximum(squares / len(rows) - means ** 2, 0)
        eligible = standard & ~excluded_mt & (detected >= 20) & (variance > 0)
        # Mean-binned dispersion rank: an explicit exploratory heuristic, not
        # a claim of verbatim Seurat/Scanpy HVG-method reproduction.
        dispersion = variance / np.maximum(means, 1e-12)
        scores = np.full(len(var), -np.inf)
        bins = pd.qcut(means[eligible], 20, labels=False, duplicates='drop')
        eligible_indices = np.flatnonzero(eligible)
        for b in np.unique(bins):
            idx = eligible_indices[np.asarray(bins) == b]
            values = np.log1p(dispersion[idx])
            scores[idx] = (values - values.mean()) / max(values.std(), 1e-12)
        selected = np.argsort(scores)[::-1][:min(hvg_count, int(eligible.sum()))]
        selected = np.sort(selected)
        if len(selected) < 3:
            raise ValueError('Too few eligible standard features')
        var = var.assign(log_cp10k_mean=means, log_cp10k_variance=variance,
                         qc_cells_detected=detected.astype(int), diagnostic_hvg=False)
        var.iloc[selected, var.columns.get_loc('diagnostic_hvg')] = True
        var.to_csv(output / 'feature_diagnostics.csv')
        state('building_sparse_hvg_matrix', cells=len(rows), features=len(selected))
        parts = []
        with h5py.File(source, 'r') as handle:
            for picked, x in blocks(handle, rows, len(var)):
                parts.append(normalized_log(x, total[picked])[:, selected])
        x = sparse.vstack(parts, format='csr')
        del parts
        sparse.save_npz(output / 'diagnostic_log_hvg.npz', x)
        atomic_json(output / 'pca_checkpoint.json', {
            'input_sha256': identity, 'cell_indices': rows.tolist(),
            'feature_indices': selected.tolist(),
            'matrix_sha256': sha256(output / 'diagnostic_log_hvg.npz'),
            'shape': list(x.shape), 'nnz': x.nnz,
            'meaning': 'diagnostic logCP10k only; not raw counts or final QC'})
        required = 256 * 1024 ** 2 + x.data.nbytes * 4
        if psutil.virtual_memory().available < required:
            raise MemoryError(f'PCA preflight needs {required} free bytes')
        state('centered_sparse_pca', cells=len(rows), features=len(selected), nnz=x.nnz)
        with threadpool_limits(limits=1):
            coordinates, loadings, eigenvalues = centered_pca(x)
        pd.DataFrame(coordinates, index=obs.index[rows],
                     columns=[f'PC{i + 1}' for i in range(coordinates.shape[1])]).to_csv(output / 'pca_coordinates.csv')
        pd.DataFrame(loadings.T, index=var.index[selected],
                     columns=[f'PC{i + 1}' for i in range(coordinates.shape[1])]).to_csv(output / 'pca_loadings.csv')
        pd.DataFrame({'pc': np.arange(len(eigenvalues)) + 1, 'variance': eigenvalues}).to_csv(output / 'pca_variance.csv', index=False)
        mapping = pd.DataFrame({'target': sorted(set(obs.loc[obs['main_analysis_stage1'], 'perturbed_gene'].astype(str)))})
        symbols = set(var.loc[standard, 'gene_name'].astype(str))
        mapping['exact_standard_gene_symbol_match'] = mapping['target'].isin(symbols)
        mapping['is_control'] = mapping['target'].eq('ctrl')
        mapping.to_csv(output / 'target_mapping_diagnostics.csv', index=False)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        diagnostic_obs = obs.iloc[rows]
        for ax, (column, title) in zip(axes, [('is_control', 'Non-targeting control'),
                                             ('total_counts', 'RNA UMI count'),
                                             ('pct_counts_mt', 'Mitochondrial RNA %')]):
            points = ax.scatter(coordinates[:, 0], coordinates[:, 1], s=2,
                                c=diagnostic_obs[column].to_numpy(dtype=float), rasterized=True)
            ax.set(xlabel='PC1', ylabel='PC2', title=title)
            fig.colorbar(points, ax=ax)
        fig.tight_layout()
        fig.savefig(output / 'expression_pca_diagnostics.png', dpi=160)
        plt.close(fig)
        outputs = {p.name: sha256(p) for p in output.iterdir() if p.suffix in {'.csv', '.png'}}
        state('complete_diagnostics_only', input_sha256=identity, cells=len(rows),
              features=len(selected), code_commit=subprocess.check_output(
                  ['git', 'rev-parse', 'HEAD'], text=True).strip(),
              doublet_assessment='pending', main_analysis_set='not_finalized',
              noncontrol_targets=len(mapping) - int(mapping['is_control'].sum()),
              targets_without_standard_symbol_match=int((~mapping['exact_standard_gene_symbol_match'] & ~mapping['is_control']).sum()),
              hvg_method='mean-binned log-CP10k dispersion heuristic; standard ENSG features; min20 cells; exclude MT',
              output_sha256=outputs)
    except Exception as exc:
        state('failed', error=repr(exc))
        raise
    finally:
        lock.unlink()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--hvg-count', type=int, default=1500)
    args = parser.parse_args()
    run(args.source, args.output, args.hvg_count)
