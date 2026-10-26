"""Export exact original plot curves from accepted derived result members.

No count matrix read or response re-estimation. Source-specific histogram bins
are preserved; a renderer may share display limits without changing the data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from perturbseq_core_io import atomic_json, identity, sha256
from perturbseq_core_effects import TYPES


def curves(records, matched, metrics, support):
    if matched.shape[0] != len(records) or not records.guide.is_unique:
        raise ValueError('Guide/response row identities differ')
    if set(metrics.guide) != set(records.guide) or not metrics.cells.eq(10).all():
        raise ValueError('Original fixed-ten-cell comparison required')
    if not metrics.groupby('guide').size().eq(8).all():
        raise ValueError('Original eight repeats required')
    if not metrics.groupby('guide').repeat.nunique().eq(8).all():
        raise ValueError('Repeated draw identity duplicated')
    finite = matched[np.isfinite(matched)]
    if not len(finite):
        raise ValueError('No original RNA effect values')
    bound = max(.05, float(np.abs(finite).max()))
    edges = np.linspace(-bound, bound, 181)
    medians = metrics.groupby('guide', sort=True).rms_log2fc.median()
    result = dict(response_features=matched.shape[1], density_edges=edges.tolist(), roles={},
                  target_coverage=[int(support.eligible_guides.eq(i).sum()) for i in range(3)]+
                                  [int(support.eligible_guides.ge(3).sum())], agreement={})
    for kind in TYPES:
        use = records.source_control_type.eq(kind).to_numpy()
        values = matched[use]
        if not len(values):
            continue
        valid = np.isfinite(values)
        weights = np.where(valid, 1/np.maximum(valid.sum(1), 1)[:, None], 0)
        hist, _ = np.histogram(values[valid], bins=edges, weights=weights[valid])
        density = hist/(hist.sum()*np.diff(edges))
        rms = np.sort(medians.loc[records.loc[use, 'guide']].to_numpy())
        if len(rms) > 1:
            x, y = np.r_[rms[0], rms], np.r_[0, np.arange(1, len(rms)+1)/len(rms)]
        else:
            x, y = rms, np.array([1.])
        result['roles'][kind] = dict(guides=int(use.sum()), density=density.tolist(),
                                    rms_x=x.tolist(), rms_y=y.tolist(),
                                    rms_median=float(np.median(rms)), rms_singleton=len(rms)==1)
    for kind in ['positive_control', 'targeting']:
        part = support.loc[support.source_control_type.eq(kind) &
                           support.assessment.eq('descriptive_assessable')]
        if len(part):
            result['agreement'][kind] = dict(targets=part.target.astype(str).tolist(),
                largest_fraction=part.largest_guide_fraction.to_numpy(float).tolist(),
                median_pair_r=part.median_pair_r.to_numpy(float).tolist(),
                guide_pairs=part.pairs.to_numpy(int).tolist())
    return result


def export(core, output, label, dataset_id):
    core, output = Path(core), Path(output)
    if (core/'owner.lock').exists():
        raise ValueError('Producer still owns this output')
    summary = json.loads((core/'summary.json').read_text())
    expected = json.loads((core/'output_manifest.json').read_text())
    expected = expected.get('outputs', expected)
    paths = {key: core/'tables'/name for key, name in [
        ('records', 'guide_effect_metadata.csv'), ('metrics', 'matched_effect_metrics.csv'),
        ('support', 'target_support.csv'), ('profiles', 'rna_effect_profiles.npz')]}
    before = {k: identity(p) for k, p in paths.items()}
    receipts = {}
    for key, path in paths.items():
        entry = expected[path.relative_to(core).as_posix()]
        if path.stat().st_size != entry['bytes']:
            raise ValueError('Accepted derived input size changed')
        if key != 'profiles' and sha256(path) != entry['sha256']:
            raise ValueError('Accepted source table hash changed')
        receipts[key] = dict(path=str(path), accepted_whole_file_sha256=entry['sha256'],
                             whole_file_rehashed_now=key!='profiles')
    records = pd.read_csv(paths['records'])
    metrics = pd.read_csv(paths['metrics'])
    support = pd.read_csv(paths['support'])
    with np.load(paths['profiles'], allow_pickle=False) as archive:
        matched = archive['matched_first']
        if not np.array_equal(archive['guides'].astype(str), records.guide.astype(str)):
            raise ValueError('Stored matched effects guide axis differs')
        genes = archive['gene_ids'].astype(str)
        if len(genes) != matched.shape[1] or not pd.Index(genes).is_unique:
            raise ValueError('Stored RNA effect feature identities differ')
    payload = curves(records, matched, metrics, support)
    if any(identity(paths[k]) != ident for k, ident in before.items()):
        raise ValueError('Derived inputs changed while exporting')
    payload.update(dataset_id=dataset_id, label=label, source_core=str(core),
        producer_code=summary.get('code', summary.get('code_commit')),
        analysis_cells=summary.get('analysis_cells', summary.get('additional_RNA_QC_cells',
                          summary.get('source_identity', {}).get('selected_cells'))),
        matched_cells=10, matched_repeats=8, guide_weighting='each_guide_equal_in_density',
        agreement_unit='one_point_per_target_median_of_guide_pairs',
        source_receipts=receipts, source_count_payload_read=False,
        consumed_matched_member_sha256=hashlib.sha256(matched.tobytes()).hexdigest(),
        consumed_gene_axis_sha256=hashlib.sha256(json.dumps(genes.tolist()).encode()).hexdigest(),
        reference_cases=support.loc[support.target.isin(['RNF7','POU5F1'])].where(pd.notna(support), None).to_dict('records'))
    # DataFrame float nulls need normal JSON nulls, not NaN.
    payload['reference_cases'] = json.loads(support.loc[support.target.isin(['RNF7','POU5F1'])].to_json(orient='records'))
    output.parent.mkdir(exist_ok=True, parents=True)
    if output.exists():
        raise FileExistsError('Do not overwrite an accepted plot export')
    atomic_json(output, payload)
    return dict(output=str(output), bytes=output.stat().st_size, sha256=sha256(output),
                source_count_payload_read=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--dataset-id', required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.core, args.output, args.label, args.dataset_id)))
