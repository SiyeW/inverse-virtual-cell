"""Tested primitives for exploratory perturbation analysis, not a full workflow.

Energy uses the empirical V-statistic from Replogle et al. 2022 (including
zero diagonals). Mann-Whitney uses tie-corrected asymptotic two-sided tests.
Cell-level p-values do not create independent biological replicates.
"""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu, false_discovery_control


def validate_pair(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1] or min(len(x), len(y), x.shape[1]) < 1:
        raise ValueError('Expected nonempty cells-by-features arrays with matching features')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite input')
    return x, y


def energy_v_statistic(x, y):
    x, y = validate_pair(x, y)
    value = 2 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()
    return max(0.0, float(value))


def energy_permutation(x, y, permutations=1000, seed=0, max_pairwise_mib=64):
    x, y = validate_pair(x, y)
    if not isinstance(permutations, int) or permutations < 1:
        raise ValueError('permutations must be a positive integer')
    n, total = len(x), len(x) + len(y)
    if total * total * 8 > max_pairwise_mib * 1024 ** 2:
        raise MemoryError('Pairwise distance budget exceeded; explicitly change sampling policy')
    combined = np.vstack([x, y])
    distances = cdist(combined, combined)
    def score(order):
        a, b = order[:n], order[n:]
        return max(0.0, float(2 * distances[np.ix_(a, b)].mean()
                        - distances[np.ix_(a, a)].mean() - distances[np.ix_(b, b)].mean()))
    observed = score(np.arange(total))
    rng = np.random.default_rng(seed)
    exceed = sum(score(rng.permutation(total)) >= observed - 1e-12
                 for _ in range(permutations))
    return observed, (exceed + 1) / (permutations + 1)


def differential_expression_block(x, control, min_cells=10, pseudocount=1.0):
    """Input: nonnegative library-size-normalized counts, not log counts.

    Returns unadjusted p-values and descriptive effects. BH must be applied
    once across ALL tested genes per contrast, not separately in each block.
    """
    x, control = validate_pair(x, control)
    if min(len(x), len(control)) < min_cells:
        raise ValueError('Insufficient cells for this exploratory contrast')
    if (x < 0).any() or (control < 0).any() or pseudocount <= 0:
        raise ValueError('Expected nonnegative normalized counts and positive pseudocount')
    statistic, p = mannwhitneyu(x, control, axis=0, alternative='two-sided', method='asymptotic')
    # A feature constant in both groups carries no contrast information.
    # Some SciPy versions emit NaN for its zero tie-corrected variance.
    constant = (np.maximum(x.max(0), control.max(0)) ==
                np.minimum(x.min(0), control.min(0)))
    statistic, p = np.asarray(statistic).copy(), np.asarray(p).copy()
    statistic[constant], p[constant] = len(x) * len(control) / 2, 1.0
    if not np.isfinite(p).all() or not np.isfinite(statistic).all():
        raise ValueError('Nonfinite test result for a nonconstant feature')
    mean_x, mean_control = x.mean(0), control.mean(0)
    return {'pvalue': p, 'rank_biserial': 2 * statistic / (len(x) * len(control)) - 1,
            'mean_normalized': mean_x, 'control_mean_normalized': mean_control,
            'log2_fold_change': np.log2((mean_x + pseudocount) / (mean_control + pseudocount)),
            'fraction_detected': (x > 0).mean(0), 'control_fraction_detected': (control > 0).mean(0)}


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    if pvalues.ndim != 1 or not len(pvalues) or not np.isfinite(pvalues).all():
        raise ValueError('Expected finite nonempty pvalue vector')
    if (pvalues < 0).any() or (pvalues > 1).any():
        raise ValueError('pvalues outside [0,1]')
    return false_discovery_control(pvalues, method='bh')
