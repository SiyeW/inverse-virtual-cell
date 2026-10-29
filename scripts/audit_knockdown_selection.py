"""Bounded feasibility audit of guide/cell RNA-repression gates on saved panels.

This is not the complete author's pipeline: the original pooled laboratory NTC
is retained, no NTC cleaning, energy test, DE test or new response profile is run.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def screen(values, reference, minimum_exclusive=25):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Invalid normalized target counts')
    if not len(values) or not np.isfinite(reference) or reference <= 0:
        return dict(reference_defined=False, guide_repression_pass=False,
                    after_cell_gate=0, final_pass=False,
                    ratio_before=None, ratio_after=None)
    before = float(values.mean() / reference)
    # Match the saved author's implementation: <0 knockdown is excluded;
    # exact equality is kept. Guide admission uses the pre-cell-filter mean.
    kept = values[values <= reference]
    guide_pass = before < 0.7
    return dict(reference_defined=True, guide_repression_pass=guide_pass,
                after_cell_gate=len(kept),
                final_pass=bool(guide_pass and len(kept) > minimum_exclusive),
                ratio_before=before,
                ratio_after=float(kept.mean()/reference) if len(kept) else None)


def summarize(frame, mask):
    x = frame.loc[mask]
    n = x.groupby('target_gene_id').size()
    return dict(guides=len(x), targets=int(x.target_gene_id.nunique()),
                targets_with_two_guides=int((n >= 2).sum()),
                original_cells=int(x.cells.sum()),
                selected_cells=int(x.after_cell_gate.sum()),
                zero_target_detected_guides=int(x.detected_cells.eq(0).sum()))


def run(config):
    output = Path(config['output'])
    output.mkdir(exist_ok=False)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if commit != config['expected_commit']:
        raise ValueError('Execution commit changed')
    summaries = []
    for d in config['datasets']:
        root = Path(d['target_result'])
        identities = dict(d['cache_metadata_sha256'])
        identities.update({'target_counts.npz': d['target_counts_sha256'],
                           'target_counts.json': d['target_counts_receipt_sha256'],
                           'guide_target_RNA.csv': d['point_table_sha256']})
        for filename, expected in identities.items():
            if sha(root/filename) != expected:
                raise ValueError('Saved input changed: '+filename)
        cells = pd.read_csv(root/'cells.csv.gz', keep_default_na=False)
        genes = pd.read_csv(root/'genes.csv.gz', keep_default_na=False)
        receipt = json.loads((root/'target_counts.json').read_text())
        if receipt['sha256'] != d['target_counts_sha256']:
            raise ValueError('Count receipt does not bind panel')
        panel = sparse.load_npz(root/'target_counts.npz').tocsc()
        positions = np.asarray(receipt['positions'], dtype=int)
        if panel.shape != (len(cells), len(positions)) or len(cells) != d['expected_cells']:
            raise ValueError('Saved axes differ')
        ntc = cells.role.eq('non_targeting_control').to_numpy()
        if int(ntc.sum()) != d['expected_NTC_cells'] or not cells.cell_barcode.is_unique:
            raise ValueError('Frozen cell cohort differs')
        if (cells.total_counts <= 0).any():
            raise ValueError('Invalid full-cell count totals')
        scale = 10000.0/cells.total_counts.to_numpy(float)
        columns = {str(genes.iloc[j].gene_id): k for k, j in enumerate(positions)}
        groups = cells.groupby('guide_id', sort=False).indices
        points = pd.read_csv(root/'guide_target_RNA.csv')
        rows = []
        maximum_error = 0.0
        for gene, sub in points.groupby('target_gene_id', sort=False):
            values = panel[:, columns[str(gene)]].toarray().ravel()*scale
            reference = float(values[ntc].mean())
            for row in sub.itertuples(index=False):
                v = values[groups.get(row.guide, np.array([], dtype=int))]
                if len(v) != int(row.cells) or int((v > 0).sum()) != int(row.detected_cells):
                    raise ValueError('Original counts or guide linkage changed')
                if len(v):
                    error = abs(float(v.mean())-row.mean_per_10000)
                    maximum_error = max(maximum_error, error)
                    if not np.isclose(v.mean(), row.mean_per_10000, rtol=1e-11, atol=1e-11):
                        raise ValueError('Original guide mean did not bridge')
                if not np.isclose(reference, row.NTC_mean_per_10000, rtol=1e-11, atol=1e-11):
                    raise ValueError('Original pooled NTC did not bridge')
                result = screen(v, reference)
                rows.append(dict(guide=row.guide, target_gene_id=gene, target=row.target,
                                 cells=int(row.cells), detected_cells=int(row.detected_cells),
                                 NTC_detected_fraction=row.NTC_detected_fraction, **result))
        frame = pd.DataFrame(rows)
        frame.to_csv(output/(d['id']+'_guide_selection.csv'), index=False)
        final = frame.final_pass
        summaries.append(dict(dataset=d['id'],
            before_count_gate=summarize(frame, frame.cells.gt(25)&frame.reference_defined),
            after_guide_and_cell_gates=summarize(frame, final),
            secondary_original_plot_detection_ge10=summarize(frame, final&frame.NTC_detected_fraction.ge(.1)),
            original_mean_bridge_max_abs_error=maximum_error,
            sources=[dict(path=str(root/n), sha256=h, bytes=(root/n).stat().st_size) for n,h in identities.items()]))
    result = dict(commit=commit, datasets=summaries,
                  all_raw_matrices_reopened=False, paper_full_pipeline=False,
                  NTC_reference='original pooled laboratory reference, unchanged',
                  not_done=['author NTC cleaning', 'batch-resolved company reference',
                            'energy and differential-expression tests', 'new full-RNA profiles'],
                  detection_ge10_is_secondary_not_author_rule=True)
    (output/'summary.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    run(json.loads(Path(parser.parse_args().config).read_text()))
