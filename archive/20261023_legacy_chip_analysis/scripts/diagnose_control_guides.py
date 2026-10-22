"""Exploratory NTC consistency in diagnostic PCA; never automatic filtering."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from analyze_expression_structure import atomic_json, sha256


def centroid_null_test(group, others, draws=1000, seed=0):
    """Finite-population label randomization, not a biological replicate test."""
    combined = np.vstack([group, others])
    n = len(group)
    observed = np.linalg.norm(group.mean(0) - others.mean(0))
    rng = np.random.default_rng(seed)
    grand_sum = combined.sum(0)
    exceed = 0
    for _ in range(draws):
        sample_sum = combined[rng.choice(len(combined), n, replace=False)].sum(0)
        value = np.linalg.norm(sample_sum / n - (grand_sum - sample_sum) / len(others))
        exceed += value >= observed
    return observed, (exceed + 1) / (draws + 1)


def run(source, diagnostics, output):
    if output.exists():
        raise FileExistsError('Do not overwrite existing control review')
    output.mkdir(parents=True)
    s = json.loads((diagnostics / 'expression_diagnostics_status.json').read_text(encoding='utf8'))
    pc_path = diagnostics / 'pca_coordinates.csv'
    if sha256(pc_path) != s['output_sha256'][pc_path.name]:
        raise ValueError('PCA hash mismatch')
    if sha256(source) != s['input_sha256']:
        raise ValueError('Source hash mismatch')
    pcs = pd.read_csv(pc_path, index_col=0)
    data = ad.read_h5ad(source, backed='r')
    obs = data.obs.loc[pcs.index].copy()
    data.file.close()
    controls = obs['is_control'].to_numpy(dtype=bool)
    ntc = pcs.loc[controls].to_numpy()
    # A common per-PC scale set by all NTCs; screen is descriptive/exploratory.
    scale = np.maximum(ntc.std(axis=0), 1e-8)
    ntc = (ntc - ntc.mean(axis=0)) / scale
    info = obs.loc[controls]
    guides = info['guide_id'].astype(str).to_numpy()
    records = []
    for i, guide in enumerate(sorted(set(guides))):
        mask = guides == guide
        row = {'guide': guide, 'cells': int(mask.sum()),
               'median_rna_umis': float(info.loc[mask, 'total_counts'].median()),
               'median_mt_pct': float(info.loc[mask, 'pct_counts_mt'].median())}
        if mask.sum() >= 20:
            row['centroid_distance'], row['label_permutation_p'] = centroid_null_test(ntc[mask], ntc[~mask], seed=i)
        records.append(row)
    table = pd.DataFrame(records)
    test = table['label_permutation_p'].notna()
    table.loc[test, 'bh_q'] = false_discovery_control(table.loc[test, 'label_permutation_p'].to_numpy())
    table.to_csv(output / 'ntc_guide_consistency.csv', index=False)
    atomic_json(output / 'summary.json', {
        'stage': 'complete_diagnostics_only', 'control_cells': int(controls.sum()),
        'observed_control_guides': len(records), 'guides_tested_min20': int(test.sum()),
        'guides_flagged_q_lt_005': int((table['bh_q'] < .05).sum()),
        'permutations': 1000, 'space': '30 diagnostic PCs scaled by NTC standard deviation',
        'input_sha256': s['input_sha256'], 'filtering_applied': False,
        'interpretation': 'Within-chip label-consistency diagnostic; not biological replicate significance; not equivalent to paper AD-based core-NTC selection',
        'output_sha256': sha256(output / 'ntc_guide_consistency.csv')})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--diagnostics', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.source, a.diagnostics, a.output)
