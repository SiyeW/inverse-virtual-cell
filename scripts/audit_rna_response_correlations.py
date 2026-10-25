"""Read-only, bounded audit of accepted guide means; no cell matrix or hit calling."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def corr(a, b):
    keep = np.isfinite(a) & np.isfinite(b)
    x, y = np.asarray(a)[keep], np.asarray(b)[keep]
    if len(x) < 3:
        return None
    x, y = x - x.mean(), y - y.mean()
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / denominator) if denominator > 1e-12 else None


def describe(values):
    x = np.array([v for v in values if v is not None and np.isfinite(v)], float)
    return dict(n=len(x), median=float(np.median(x)) if len(x) else None,
                q25=float(np.quantile(x, .25)) if len(x) else None,
                q75=float(np.quantile(x, .75)) if len(x) else None)


class BudgetReader:
    def __init__(self, stream, budget, deadline):
        self.stream, self.budget, self.deadline, self.bytes_read = stream, budget, deadline, 0

    def read(self, size=-1):
        if time.monotonic() > self.deadline:
            raise TimeoutError('Profile read deadline reached; no automatic retry')
        if size < 0:
            position = self.stream.tell()
            size = self.stream.seek(0, 2) - position
            self.stream.seek(position)
        if size > self.budget - self.bytes_read:
            raise ValueError('Profile read budget exceeded')
        value = self.stream.read(size)
        self.bytes_read += len(value)
        if time.monotonic() > self.deadline:
            raise TimeoutError('Profile read deadline reached; no automatic retry')
        return value

    def __getattr__(self, name):
        return getattr(self.stream, name)


def identity(path):
    s = path.stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)


def matched_other(records, a, b):
    """Keep guide a, replace b by one other-target guide of the same role and nearest n."""
    target, role, count = records.iloc[b][['target', 'source_control_type', 'cells']]
    candidates = records.index[(records.target != target) &
                               (records.source_control_type == role)].tolist()
    if not candidates:
        return None
    return min(candidates, key=lambda j: (
        abs(np.log(records.iloc[j].cells / count)),
        hashlib.sha256((str(records.iloc[a].guide) + ':' + str(records.iloc[b].guide) + ':' +
                        str(records.iloc[j].guide)).encode()).hexdigest()))


def calculate(records, pairs, means, effects, guides, genes):
    if not np.array_equal(np.asarray(guides, str), records.guide.astype(str)):
        raise ValueError('Guide order mismatch')
    if means.shape != effects.shape or means.shape != (len(records), len(genes)):
        raise ValueError('Profile axes mismatch')
    if len(set(genes)) != len(genes) or not np.isfinite(means).all() or (means < 0).any():
        raise ValueError('Invalid feature identities or means')
    ntc = records.source_control_type.eq('non_targeting_control').to_numpy()
    counts = records.cells.to_numpy(float)
    total_ntc = counts[ntc].sum()
    reference = np.average(means[ntc], axis=0, weights=counts[ntc])
    expected = np.log2((means + 1) / (reference + 1))
    for i in np.flatnonzero(ntc):
        center = (reference * total_ntc - means[i] * counts[i]) / (total_ntc - counts[i])
        expected[i] = np.log2((means[i] + 1) / (center + 1))
    for i, value in enumerate(records.direct_feature_index):
        if pd.notna(value):
            expected[i, int(value)] = np.nan
    if not np.array_equal(np.isnan(expected), np.isnan(effects)):
        raise ValueError('Direct-feature mask differs')
    np.testing.assert_allclose(expected, effects, atol=1e-10, rtol=1e-10, equal_nan=True)
    lookup = {g: i for i, g in enumerate(records.guide)}
    filters = {'all_accepted': np.ones(len(genes), bool),
               **{f'ntc_mean_ge_{n}': reference >= n for n in (1, 5, 10)}}
    rows = []
    for pair in pairs.itertuples():
        a, b = lookup[pair.guide_a], lookup[pair.guide_b]
        full_r = corr(effects[a], effects[b])
        if full_r is None or abs(full_r - pair.pearson) > 1e-10:
            raise ValueError('Accepted pair correlation cannot be reproduced')
        other = matched_other(records, a, b)
        for label, keep in filters.items():
            rows.append(dict(target=pair.target, guide_a=pair.guide_a, guide_b=pair.guide_b,
                role=records.iloc[a].source_control_type, filter=label, features=int(keep.sum()),
                minimum_cells=min(pair.cells_a, pair.cells_b),
                baseline_mean_r=corr(means[a, keep], means[b, keep]),
                baseline_log1p_r=corr(np.log1p(means[a, keep]), np.log1p(means[b, keep])),
                response_r=corr(effects[a, keep], effects[b, keep]),
                different_target_guide=None if other is None else records.iloc[other].guide,
                different_target_cells=None if other is None else int(counts[other]),
                replaced_guide_cells=int(counts[b]),
                different_target_r=None if other is None else corr(effects[a, keep], effects[other, keep])))
    table = pd.DataFrame(rows)
    # Example chosen by coverage alone, never by the resulting correlation.
    best = pairs.assign(minimum=pairs[['cells_a', 'cells_b']].min(axis=1)).sort_values(
        ['minimum', 'target', 'guide_a', 'guide_b'], ascending=[False, True, True, True]).iloc[0]
    a, b = lookup[best.guide_a], lookup[best.guide_b]
    example = pd.DataFrame(dict(gene_id=genes, ntc_mean=reference, mean_a=means[a], mean_b=means[b],
                                log2fc_a=effects[a], log2fc_b=effects[b]))
    summary = dict(pairs_reproduced=len(pairs), eligible_guides=len(records), features=len(genes),
        ntc_cells=int(total_ntc), ntc_guides=int(ntc.sum()),
        response_formula='log2((mean depth-normalized RNA + 1)/(NTC mean + 1)); whole-guide leave-out for NTC',
        feature_filters_predeclared=['accepted feature set', 'NTC mean >= 1', 'NTC mean >= 5', 'NTC mean >= 10'],
        different_target_background='one deterministic same-role nearest-cell-count replacement for guide b; a retained; descriptive background, not a randomized null',
        example_selection='largest minimum cell count, target/guide lexical tie-break; no correlation selection',
        example=dict(target=str(best.target), guide_a=str(best.guide_a), guide_b=str(best.guide_b),
                     cells_a=int(best.cells_a), cells_b=int(best.cells_b)),
        by_filter={label:dict(features=int(part.features.iloc[0]),
            baseline_mean=describe(part.baseline_mean_r), baseline_log1p=describe(part.baseline_log1p_r),
            response=describe(part.response_r), different_target=describe(part.different_target_r),
            paired_response_minus_different=describe(part.response_r - part.different_target_r))
            for label, part in table.groupby('filter', sort=False)},
        by_minimum_cells={str(n):describe(table.loc[table['filter'].eq('all_accepted') &
            table.minimum_cells.ge(n), 'response_r']) for n in (10, 20, 50, 100, 200)},
        no_functional_hit_calling=True, original_filters_changed=False, per_cell_matrix_read=False,
        feature_selection_uses_guide_contrasts=False, independent_biological_replicates=False)
    return table, example, summary


def run(core, output, manifest_sha, code_commit, read_mib=80, seconds=90):
    started = time.monotonic()
    core, output = Path(core).resolve(), Path(output).resolve()
    if output.exists() or output == core or output.is_relative_to(core) or core.is_relative_to(output):
        raise ValueError('Distinct unused audit output required')
    if (core/'owner.lock').exists():
        raise ValueError('Accepted core has an active writer')
    raw = (core/'output_manifest.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest_sha:
        raise ValueError('Core manifest differs')
    manifest = json.loads(raw)
    for name in ['tables/guide_effect_metadata.csv', 'tables/within_target_pairs.csv']:
        if sha(core/name) != manifest[name]['sha256']:
            raise ValueError('Accepted small table differs')
    records = pd.read_csv(core/'tables/guide_effect_metadata.csv')
    pairs = pd.read_csv(core/'tables/within_target_pairs.csv')
    path = core/'tables/rna_effect_profiles.npz'
    before = identity(path)
    if before['bytes'] != manifest['tables/rna_effect_profiles.npz']['bytes']:
        raise ValueError('Profile file size differs')
    with path.open('rb') as stream:
        reader = BudgetReader(stream, read_mib * 1024**2, time.monotonic() + seconds)
        with np.load(reader, allow_pickle=False) as values:
            guides, genes = values['guides'].astype(str), values['gene_ids'].astype(str)
            means, effects = values['means'], values['log2fc']
    after = identity(path)
    if any(before[k] != after[k] for k in ['bytes', 'mtime_ns', 'device']):
        raise ValueError('Accepted profile changed during read')
    table, example, summary = calculate(records, pairs, means, effects, guides, genes)
    summary.update(stage='complete', code_commit=code_commit, helper_sha256=sha(__file__),
        core_manifest_sha256=manifest_sha, source_npz_full_sha256_reused=manifest['tables/rna_effect_profiles.npz']['sha256'],
        source_npz_fresh_full_hash=False, selected_zip_members_crc_checked=True,
        source_identity_before=before, source_identity_after=after, compressed_profile_bytes_read=reader.bytes_read,
        read_budget_bytes=read_mib*1024**2, elapsed_seconds=time.monotonic()-started)
    output.mkdir()
    table.to_csv(output/'pair_comparisons.csv', index=False)
    example.to_csv(output/'example_profiles.csv', index=False)
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf8')
    proofs = {p.name:dict(bytes=p.stat().st_size, sha256=sha(p)) for p in output.iterdir()}
    (output/'output_manifest.json').write_text(json.dumps(proofs, indent=2), encoding='utf8')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--core', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--manifest-sha256', required=True)
    p.add_argument('--code-commit', required=True)
    p.add_argument('--read-mib', type=int, default=80)
    p.add_argument('--read-seconds', type=int, default=90)
    a = p.parse_args()
    with threadpool_limits(limits=2):
        print(json.dumps(run(a.core, a.output, a.manifest_sha256, a.code_commit, a.read_mib, a.read_seconds), indent=2))
