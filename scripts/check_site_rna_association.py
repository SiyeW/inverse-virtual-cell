"""Exploratory RNA-state association of observed RNA indels, never a KO caller.

Energy U statistics use the frozen normalized 2000-RNA panel, excluding the
target RNA. Evidence labels are permuted within sgRNA, preserving guide-specific
label frequencies. RNA-amount residualization and >=2-molecule support are
reported as sensitivity analyses; no cells or sgRNAs are selected by the tests.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


def energy_from_indicators(distance, labels):
    labels = np.atleast_2d(labels).astype(float)
    n = distance.shape[0]
    na = labels.sum(axis=1)
    nb = n - na
    if np.any(na < 2) or np.any(nb < 2):
        raise ValueError('Energy U requires at least two observations per group')
    within_a = np.einsum('ij,ij->i', labels @ distance, labels)
    rows_a = labels @ distance.sum(axis=1)
    between = rows_a - within_a
    within_b = distance.sum() - 2 * rows_a + within_a
    return 2 * between / (na * nb) - within_a / (na * (na - 1)) - within_b / (nb * (nb - 1))


def within_guide_permutations(labels, guides, repeats=999, seed=20261103):
    rng = np.random.default_rng(seed)
    result = np.repeat(np.asarray(labels)[None, :], repeats, axis=0)
    mutable = 0
    for g in np.unique(guides):
        idx = np.flatnonzero(guides == g)
        if len(np.unique(labels[idx])) < 2:
            continue
        mutable += len(idx)
        for i in range(repeats):
            result[i, idx] = rng.permutation(labels[idx])
    return result, mutable


def adjusted_bh(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    out = np.empty(len(values))
    out[order] = np.minimum(1, np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    return out


def run(cfg):
    import anndata as ad
    from threadpoolctl import threadpool_limits
    import os, socket
    assert os.getuid() == 2601 and socket.gethostname() == 'rhesusbase'
    if hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
    sequence_root = Path(cfg['sequence_output'])
    out = Path(cfg['output']) / 'association'
    out.mkdir(parents=True, exist_ok=True)
    selected = json.loads(Path(cfg['selected_guides']).read_text())
    targets = sorted({r['target_gene_id'] for r in selected})
    catalogue = pd.read_csv(cfg['catalogue'], keep_default_na=False).set_index('guide')
    for dataset in cfg['datasets']:
        a = ad.read_h5ad(dataset['RNA_checkpoint'])
        meta = a.obs.copy()
        meta['target_gene_id'] = meta.guide_id.astype(str).map(catalogue.target_gene_id)
        sequence = pd.DataFrame(json.loads((sequence_root / 'sequence' / (dataset['id'] + '.json')).read_text())['cells'])
        assert not sequence.duplicated(['guide', 'cell_barcode']).any()
        evidence = sequence.set_index(['guide', 'cell_barcode'])
        results = []
        for target in targets:
            ids = np.flatnonzero(meta.target_gene_id.eq(target).to_numpy())
            cells = meta.iloc[ids]
            joined = evidence.reindex(pd.MultiIndex.from_arrays([cells.guide_id.astype(str), cells.index.astype(str)], names=['guide', 'cell_barcode'])).fillna(0)
            for kind, minimum in [('any_indel_RNA', 1), ('two_or_more_indel_RNA_molecules', 2)]:
                positive = joined.indel_families.to_numpy() >= minimum
                negative = ((joined.indel_families.to_numpy() == 0) & (joined.SNV_families.to_numpy() == 0)
                            & (joined.cut_ref_families.to_numpy() > 0) & (joined.ambiguous_families.to_numpy() == 0))
                use = positive | negative
                labels = positive[use].astype(int)
                na, nb = int(positive.sum()), int(negative.sum())
                row = dict(dataset=dataset['id'], target_gene_id=target,
                           target=catalogue.loc[catalogue.target_gene_id.eq(target), 'target_gene_name'].iloc[0],
                           comparison=kind, RNA_records=len(cells), indel_RNA_records=na,
                           reference_window_records=nb, compared_records=int(use.sum()),
                           no_comparison_records=int((~use).sum()), status='insufficient_observed_records')
                if na >= 2 and nb >= 2:
                    guides = cells.guide_id.astype(str).to_numpy()[use]
                    seed = int.from_bytes((dataset['id'] + target).encode()[-4:], 'little') + 20261103
                    permutations, mutable = within_guide_permutations(labels, guides, seed=seed)
                    x = a.X[ids[use]][:, np.asarray(a.var_names != target)].toarray().astype(float)
                    amount = np.log1p(cells.total_counts.to_numpy(dtype=float)[use])
                    z = np.column_stack([np.ones(len(amount)), amount - amount.mean()])
                    with threadpool_limits(limits=2):
                        residual = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
                        for label, matrix in [('raw', x), ('RNA_amount_residual', residual)]:
                            d = cdist(matrix, matrix)
                            observed = float(energy_from_indicators(d, labels)[0])
                            null = energy_from_indicators(d, permutations)
                            row['energy_' + label] = observed
                            row['p_' + label] = float((1 + np.count_nonzero(null >= observed - 1e-12)) / (len(null) + 1)) if mutable else np.nan
                    row.update(status='tested_record_association' if mutable else 'no_within_guide_label_exchange',
                               mutable_records=mutable, indel_guides=int(len(np.unique(guides[labels == 1]))),
                               log_RNA_difference=float(amount[labels == 1].mean() - amount[labels == 0].mean()))
                results.append(row)
            print(json.dumps(dict(stage='association', dataset=dataset['id'], target=target)), flush=True)
        result = pd.DataFrame(results)
        for kind in result.comparison.unique():
            for metric in ['raw', 'RNA_amount_residual']:
                pcol = 'p_' + metric
                if pcol not in result:
                    result[pcol] = np.nan
                valid = result.comparison.eq(kind) & result[pcol].notna()
                result.loc[valid, 'q_' + metric] = adjusted_bh(result.loc[valid, pcol])
        result.to_csv(out / (dataset['id'] + '.csv'), index=False)
    (out / 'method.json').write_text(json.dumps(dict(
        code_commit=cfg['commit'], sequence_output=cfg['sequence_output'],
        method='multivariate Euclidean energy U on full frozen log-normalized RNA panel, target RNA excluded',
        reference_state='RNA cut-window reference support with no observed indel; not DNA wild type',
        permutations=999, permutation_blocks='sgRNA; counts preserved within each guide',
        multiple_comparisons='BH separately within dataset and molecular-support sensitivity',
        RNA_amount_sensitivity='project each RNA off intercept and log1p(authoritative total RNA), without using evidence label',
        purpose='exploratory association only; no KO calling or downstream cell selection',
        completed_at=time.strftime('%Y-%m-%dT%H:%M:%S%z')), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--config', type=Path, required=True)
    run(json.loads(p.parse_args().config.read_text()))
