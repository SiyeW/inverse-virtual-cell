"""Diagnostic Scrublet scores, not an automatic cell-exclusion decision."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from sklearn import config_context
from threadpoolctl import threadpool_limits

from analyze_expression_structure import atomic_json, blocks, sha256


def run(source, checkpoint_dir, output, limit=0, mean_center=False, sim_ratio=1.0):
    output.mkdir(parents=True, exist_ok=True)
    status = output / 'status.json'
    if status.exists():
        raise FileExistsError('Inspect existing run before restarting')
    lock = output / 'owner.lock'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    started = time.monotonic()
    def state(stage, **extra):
        atomic_json(status, {'stage': stage, 'pid': os.getpid(),
                    'elapsed_seconds': time.monotonic() - started, **extra})
    try:
        state('verify_input')
        checkpoint = json.loads((checkpoint_dir / 'pca_checkpoint.json').read_text(encoding='utf8'))
        if sha256(source) != checkpoint['input_sha256']:
            raise ValueError('Input counts hash mismatch')
        rows = np.asarray(checkpoint['cell_indices'], dtype=int)
        columns = np.asarray(checkpoint['feature_indices'], dtype=int)
        if limit:
            rows = np.sort(np.random.default_rng(0).choice(rows, min(limit, len(rows)), replace=False))
        required_mib = 700 if limit else (2000 if mean_center else 1400)
        if psutil.virtual_memory().available < required_mib * 1024**2:
            raise MemoryError('Insufficient free memory for diagnostic Scrublet')
        a = ad.read_h5ad(source, backed='r')
        obs, var = a.obs.iloc[rows].copy(), a.var.iloc[columns].copy()
        n_vars = a.n_vars
        a.file.close()
        state('building_raw_count_subset', cells=len(rows), features=len(columns))
        parts = []
        with h5py.File(source, 'r') as handle:
            for _, counts in blocks(handle, rows, n_vars):
                parts.append(counts[:, columns])
        counts = sparse.vstack(parts, format='csr')
        del parts
        data = ad.AnnData(counts, obs=obs, var=var)
        import scanpy as sc
        state('scrublet', cells=len(rows), features=len(columns))
        # The uncentered 1:1 option bounds memory; centered 2:1 sensitivity
        # uses Scrublet's default embedding/simulation choices, when safe.
        with config_context(working_memory=64), threadpool_limits(limits=1):
            sc.pp.scrublet(data, expected_doublet_rate=0.05, sim_doublet_ratio=sim_ratio,
                          mean_center=mean_center, n_prin_comps=30,
                          use_approx_neighbors=False, random_state=0)
        data.obs.to_csv(output / 'cell_doublet_diagnostics.csv.gz')
        result = data.uns['scrublet']
        simulated = np.asarray(result['doublet_scores_sim'])
        observed = data.obs['doublet_score'].to_numpy()
        if not np.isfinite(observed).all() or not np.isfinite(simulated).all():
            raise ValueError('Nonfinite scores')
        np.save(output / 'simulated_doublet_scores.npy', simulated)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(observed, bins=60, alpha=.6, density=True, label='Observed cells')
        ax.hist(simulated, bins=60, alpha=.6, density=True, label='Simulated doublets')
        threshold = result.get('threshold')
        if threshold is not None:
            ax.axvline(threshold, color='black', linestyle='--', label='Automatic candidate threshold')
        ax.set(xlabel='Diagnostic doublet score', ylabel='Density')
        ax.legend()
        fig.tight_layout()
        fig.savefig(output / 'doublet_score_distribution.png', dpi=160)
        plt.close(fig)
        state('complete_smoke_only' if limit else 'complete_diagnostics_only',
              cells=len(rows), input_sha256=checkpoint['input_sha256'],
              code_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
              expected_rate_assumption=0.05, experimental_rate='unverified',
              sim_doublet_ratio=sim_ratio, mean_center=mean_center, neighbors='exact; working_memory64MiB',
              candidate_threshold=float(threshold) if threshold is not None else None,
              candidate_calls=int(data.obs['predicted_doublet'].sum()) if threshold is not None else None,
              filtering_applied=False, source_counts_modified=False,
              output_sha256={p.name:sha256(p) for p in output.iterdir()
                             if p.suffix in {'.gz','.npy','.png'}})
    except Exception as exc:
        state('failed', error=repr(exc))
        raise
    finally:
        lock.unlink()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--checkpoint-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--mean-center',action='store_true')
    p.add_argument('--sim-ratio',type=float,default=1.0)
    args = p.parse_args()
    run(args.source,args.checkpoint_dir,args.output,args.limit,args.mean_center,args.sim_ratio)
