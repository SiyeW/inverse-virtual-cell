"""Bounded, counts-preserving supplemental QC primitives.

These diagnostics do not select cells, infer DNA editing, replace the DE
method, or turn sgRNA pairs into biological replicates.
"""
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, rankdata


def distribution_summary(values, thresholds=()):
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or not len(x) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('Expected nonempty finite nonnegative observations')
    q10, med, q90 = np.quantile(x, [.1, .5, .9], method='linear')
    mean = float(x.mean())
    ordered = np.sort(x)
    gini = (2*np.dot(np.arange(1, len(x)+1), ordered)/(len(x)*x.sum())
            - (len(x)+1)/len(x)) if mean > 0 else None
    result = dict(n=len(x), minimum=float(x.min()), p10=float(q10), median=float(med),
                  p90=float(q90), maximum=float(x.max()), mean=mean,
                  mad=float(np.median(np.abs(x-med))),
                  cv_population=float(x.std(ddof=0)/mean) if mean > 0 else None,
                  gini=gini, zero_n=int((x == 0).sum()), zero_fraction=float((x == 0).mean()),
                  p90_over_p10=float(q90/q10) if q10 > 0 else None,
                  skew_status='defined' if q10 > 0 else 'undefined_zero_p10')
    for cutoff in thresholds:
        result[f'n_below_{cutoff}'] = int((x < cutoff).sum())
        result[f'fraction_below_{cutoff}'] = float((x < cutoff).mean())
    return result


def validate_guide_coverage(obs, coverage):
    required = {'name', 'target_gene_name', 'single_guide_before_qc',
                'single_guide_after_stage1_qc'}
    if not required.issubset(coverage.columns) or coverage.name.duplicated().any():
        raise ValueError('Invalid frozen guide reference/coverage table')
    frame = coverage.copy()
    frame['guide_type'] = np.where(frame.target_gene_name.eq('Non-Targeting'), 'NTC', 'targeting')
    single = obs.guide_status.astype(str).eq('single_guide')
    for field, mask in [('single_guide_before_qc', single),
                        ('single_guide_after_stage1_qc', obs.main_analysis_stage1.astype(bool))]:
        observed = obs.loc[mask, 'guide_id'].astype(str).value_counts()
        if not set(observed.index).issubset(set(frame.name.astype(str))):
            raise ValueError('Observed guide missing from library')
        expected = frame.name.astype(str).map(observed).fillna(0).to_numpy(dtype=int)
        if not np.array_equal(expected, frame[field].to_numpy(dtype=int)):
            raise ValueError('Guide coverage does not match source cell assignments')
    return frame


def plan_null_splits(n_ntc, sizes, control_cells, repeats, seed):
    if min(n_ntc, control_cells, repeats) < 1 or any(n < 1 or n+control_cells > n_ntc for n in sizes):
        raise ValueError('Insufficient distinct NTCs for disjoint pseudo-treatment/control')
    rng = np.random.default_rng(seed)
    result = []
    for size in sizes:
        for replicate in range(repeats):
            order = rng.choice(n_ntc, size+control_cells, replace=False)
            result.append(dict(comparison=f'n{size}_r{replicate+1:02d}', cells=int(size),
                               replicate=replicate+1,
                               treatment=np.sort(order[:size]), control=np.sort(order[size:])))
    return result


def library_normalize(counts, totals, reference_total):
    totals = np.asarray(totals, dtype=float)
    if counts.shape[0] != len(totals) or (totals <= 0).any() or not np.isfinite(totals).all():
        raise ValueError('Invalid source row totals')
    if not np.isfinite(reference_total) or reference_total <= 0:
        raise ValueError('Invalid fixed library-size normalization target')
    if (counts.data < 0).any():
        raise ValueError('Negative source counts')
    return (sparse.diags(reference_total/totals) @ counts.astype(float)).tocsr()


def masked_response_profiles(normalized, guides, targets, center, scale, feature_names,
                             min_cells=10):
    guides, targets = np.asarray(guides, dtype=str), np.asarray(targets, dtype=str)
    center, scale = np.asarray(center, dtype=float), np.asarray(scale, dtype=float)
    names = np.asarray(feature_names, dtype=str)
    if normalized.shape != (len(guides), len(center)) or len(targets) != len(guides):
        raise ValueError('Profile identity/shape mismatch')
    if names.shape != scale.shape or scale.shape != center.shape or not np.isfinite(scale).all() or (scale <= 0).any():
        raise ValueError('Profiles require finite positive frozen control SD')
    records, profiles = [], []
    for guide in sorted(set(guides)):
        mask = guides == guide
        labels = set(targets[mask])
        if len(labels) != 1:
            raise ValueError('One guide maps to multiple perturbation targets')
        target = next(iter(labels))
        if int(mask.sum()) < min_cells:
            continue
        mean = np.asarray(normalized[mask].mean(axis=0)).ravel()
        response = (mean-center)/scale
        direct = np.flatnonzero(names == target)
        masked = len(direct) == 1
        if masked:
            response[direct[0]] = 0
        records.append(dict(guide=guide, target=target, cells=int(mask.sum()),
                            guide_type='NTC' if target == 'ctrl' else 'targeting',
                            direct_target_feature_masked=masked))
        profiles.append(response)
    if len(profiles) < 2:
        raise ValueError('Insufficient assessable guide profiles')
    return pd.DataFrame(records), np.asarray(profiles)


def profile_correlations(profiles):
    values = np.asarray(profiles, dtype=float)
    if values.ndim != 2 or values.shape[1] < 2 or not np.isfinite(values).all():
        raise ValueError('Invalid response profiles')
    centered = values-values.mean(axis=1, keepdims=True)
    lengths = np.linalg.norm(centered, axis=1)
    valid = lengths > 1e-12
    unit = np.zeros_like(centered)
    unit[valid] = centered[valid]/lengths[valid, None]
    corr = np.clip(unit @ unit.T, -1, 1)
    corr[~valid, :] = np.nan
    corr[:, ~valid] = np.nan
    return corr, valid


def guide_agreement(records, profiles, max_between=10000, seed=0):
    corr, valid = profile_correlations(profiles)
    rng = np.random.default_rng(seed)
    within, ntc, between = [], [], []
    for i, j in combinations(range(len(records)), 2):
        if not valid[i] or not valid[j]:
            continue
        a, b = records.iloc[i], records.iloc[j]
        row = dict(guide_a=a.guide, guide_b=b.guide, target_a=a.target, target_b=b.target,
                   cells_a=int(a.cells), cells_b=int(b.cells), min_cells=int(min(a.cells, b.cells)),
                   response_correlation=float(corr[i, j]))
        if a.guide_type == b.guide_type == 'NTC':
            ntc.append(row | dict(pair_type='NTC-NTC'))
        elif a.guide_type == b.guide_type == 'targeting':
            if a.target == b.target:
                within.append(row | dict(pair_type='same-target'))
            else:
                between.append(row | dict(pair_type='different-target'))
    if len(between) > max_between:
        between = [between[i] for i in np.sort(rng.choice(len(between), max_between, replace=False))]
    columns = ['guide_a', 'guide_b', 'target_a', 'target_b', 'cells_a', 'cells_b',
               'min_cells', 'response_correlation', 'pair_type']
    return pd.DataFrame(within+ntc+between, columns=columns), valid


def null_test_summary(pvalues, minimum=1e-300):
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not len(p) or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Invalid full testing-family p-values')
    q = false_discovery_control(p, method='bh')
    return dict(tested_features=len(p), raw_p_lt_005=int((p < .05).sum()),
                raw_p_lt_001=int((p < .01).sum()), bh_q_lt_005=int((q < .05).sum()),
                min_p=float(p.min()), min_q=float(q.min()),
                median_neg_log10_p=float(np.median(-np.log10(np.maximum(p, minimum)))))


def permutation_rank_diagnostic(x, y, permutations=19999, seed=0):
    """Bounded independent rank-sum randomization check, not a DE replacement.

    Uses a declared two-sided absolute-deviation statistic conditional on all
    observed values, preserving ties and unequal sample sizes. Exact ordering
    and count data are unchanged; Monte Carlo p uses the plus-one correction.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or not min(len(x), len(y)) or permutations < 1:
        raise ValueError('Nonempty one-dimensional samples and positive budget required')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite diagnostic observations')
    if len(x) > len(y): x, y = y, x
    n, m = len(x), len(y)
    ranks = rankdata(np.concatenate([x, y]), method='average')
    observed_u = ranks[:n].sum()-n*(n+1)/2
    distance = abs(observed_u-n*m/2)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        selected = rng.choice(n+m, n, replace=False)
        u = ranks[selected].sum()-n*(n+1)/2
        exceed += abs(u-n*m/2) >= distance-1e-12
    return dict(pvalue=(int(exceed)+1)/(permutations+1), exceedances=int(exceed),
                permutations=permutations, statistic=float(distance),
                statistic_definition='abs(U - n*m/2)', seed=seed,
                interpretation='A spot diagnostic only; not recalibrated genome-wide DE')
