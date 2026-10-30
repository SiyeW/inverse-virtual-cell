"""Full-population energy statistics and independent/restricted permutations.

The V statistic matches released pertpy Edistance. The U effect is the usual
unequal-size unbiased estimate used in the accepted exploration. They are not
assumed to have identical permutation ranks at unequal group sizes.
"""
import numpy as np
from numba import njit, prange


def energy_from_sums(n, m, inside, selected_rows, total):
    if min(n, m) < 2:
        raise ValueError('Both groups need at least two observations')
    cross = selected_rows - inside
    other = total + inside - 2 * selected_rows
    v = 2 * cross / (n*m) - inside/(n*n) - other/(m*m)
    u = 2 * cross / (n*m) - inside/(n*(n-1)) - other/(m*(m-1))
    return float(v), float(u)


@njit(cache=True, parallel=True)
def subset_statistics(samples, control, cross, query, row_sums, total):
    """Exact selected-pair sums; no large group-mask matrix is allocated."""
    nc = control.shape[0]
    n = samples.shape[1]
    m = row_sums.size - n
    output = np.empty((samples.shape[0], 2), dtype=np.float64)
    for r in prange(samples.shape[0]):
        within = 0.0
        selected_rows = 0.0
        for i in range(n):
            a = samples[r, i]
            selected_rows += row_sums[a]
            for j in range(i):
                b = samples[r, j]
                if a < nc and b < nc:
                    value = control[a, b]
                elif a >= nc and b >= nc:
                    value = query[a-nc, b-nc]
                elif a >= nc:
                    value = cross[a-nc, b]
                else:
                    value = cross[b-nc, a]
                within += 2 * value
        between = selected_rows - within
        other = total + within - 2*selected_rows
        output[r, 0] = 2*between/(n*m)-within/(n*n)-other/(m*m)
        output[r, 1] = 2*between/(n*m)-within/(n*(n-1))-other/(m*(m-1))
    return output


def permutation_samples(size, observed, permutations, seed, strata=None):
    """Uniform fixed-size label allocations, optionally within fixed blocks."""
    observed = np.asarray(observed, dtype=np.int64)
    if len(np.unique(observed)) != len(observed) or min(len(observed), size-len(observed)) < 2:
        raise ValueError('Invalid observed allocation')
    rng = np.random.default_rng(seed)
    result = np.empty((permutations, len(observed)), dtype=np.int32)
    if strata is None:
        for i in range(permutations):
            result[i] = np.sort(rng.choice(size, len(observed), replace=False))
    else:
        strata = np.asarray(strata)
        if strata.shape != (size,):
            raise ValueError('Strata do not cover the full comparison cohort')
        blocks=[]
        for level in np.unique(strata):
            pool=np.flatnonzero(strata == level)
            count=np.count_nonzero(strata[observed] == level)
            if count: blocks.append((pool, count))
        for i in range(permutations):
            offset=0
            for pool,count in blocks:
                result[i,offset:offset+count]=rng.choice(pool,count,replace=False)
                offset+=count
            result[i].sort()
    return result


def depth_strata(total_rna, bins=10):
    """Label-blind RNA-depth quantiles of the complete compared population."""
    values=np.asarray(total_rna, dtype=np.float64)
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('Invalid RNA totals')
    cuts=np.unique(np.quantile(np.log1p(values), np.arange(1,bins)/bins))
    return np.searchsorted(cuts,np.log1p(values),side='right')


def depth_complexity_strata(total_rna, detected_rna, bins=10):
    """Joint, label-blind quantile blocks of two standard count covariates."""
    total_rna=np.asarray(total_rna);detected_rna=np.asarray(detected_rna)
    if total_rna.shape!=detected_rna.shape:raise ValueError('Covariate axes differ')
    a=depth_strata(total_rna,bins);b=depth_strata(detected_rna,bins)
    return np.unique(a*bins+b,return_inverse=True)[1]


def permutation_pvalue(observed, null):
    null=np.asarray(null)
    tolerance=1e-10*max(1.,abs(float(observed)))
    extreme=int(np.count_nonzero(null >= observed-tolerance))
    return (extreme+1)/(len(null)+1), extreme
