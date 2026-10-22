"""Numerical primitives for the recorded Replogle-2022-aligned workflow.

AD implements the two-sample midrank statistic (Scholz & Stephens 1987),
validated against scipy.stats.anderson_ksamp. Tail probabilities are supplied
by a separately provenance-tracked grid from R kSamples::ad.pval(m=1,version=1),
not SciPy's capped tabulation. No laboratory data or external R code lives here.
"""
from functools import lru_cache
from pathlib import Path
import numpy as np
from numba import njit
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import cdist, pdist
from scipy.stats import false_discovery_control


def reference_de_domain(pvalues, relative_z_defined):
    """BH on the finite positive control-z domain; store others as untested."""
    p=np.asarray(pvalues,dtype=float);mask=np.asarray(relative_z_defined,dtype=bool)
    if p.ndim!=1 or mask.shape!=p.shape or not mask.any():
        raise ValueError('Expected a nonempty defined control-z testing domain')
    if not np.isfinite(p[mask]).all() or np.any((p[mask]<0)|(p[mask]>1)):
        raise ValueError('Invalid p values in the actual testing domain')
    retained=np.full(p.shape,np.nan);adjusted=np.full(p.shape,np.nan)
    retained[mask]=p[mask];adjusted[mask]=false_discovery_control(p[mask])
    return retained,adjusted


@lru_cache(maxsize=1024)
def ad_standard_deviation(n, m):
    total = n + m
    if min(n, m) < 1 or total < 4:
        raise ValueError('AD standardized variance requires positive samples and N>=4')
    reciprocal = np.cumsum(1.0 / np.arange(total-1, 1, -1))
    h = reciprocal[-1] + 1.0
    g = np.sum(reciprocal / np.arange(2, total))
    H, k = 1/n + 1/m, 2
    a = (4*g-6)*(k-1)+(10-6*g)*H
    b = (2*g-4)*k*k+8*h*k+(2*g-14*h-4)*H-8*h+4*g-6
    c = (6*h+2*g-2)*k*k+(4*h-4*g+6)*k+(2*h-6)*H+4*h
    d = (2*h+6)*k*k-4*h*k
    variance = (a*total**3+b*total**2+c*total+d)/((total-1)*(total-2)*(total-3))
    if not np.isfinite(variance) or variance <= 0:
        raise ValueError('Invalid AD variance')
    return float(np.sqrt(variance))


@njit(cache=True)
def _ad_sorted_columns(x, y, deviation):
    """Merge sorted samples by tied runs; each run contributes its midrank."""
    n, m, genes = len(x), len(y), x.shape[1]
    total = n + m
    result = np.empty(genes)
    for column in range(genes):
        if x[0,column] == x[-1,column] == y[0,column] == y[-1,column]:
            result[column] = -np.inf  # no information, caller assigns p=1
            continue
        ix, iy, accum = 0, 0, 0.0
        while ix < n or iy < m:
            if iy == m or (ix < n and x[ix,column] < y[iy,column]):
                value = x[ix,column]
            else:
                value = y[iy,column]
            endx, endy = ix, iy
            while endx < n and x[endx,column] == value: endx += 1
            while endy < m and y[endy,column] == value: endy += 1
            freq = endx-ix + endy-iy
            rank = ix+iy+freq/2.0
            denominator = rank*(total-rank)-total*freq/4.0
            if denominator > 0:
                midx, midy = (ix+endx)/2.0, (iy+endy)/2.0
                accum += freq/total*((total*midx-n*rank)**2/n + (total*midy-m*rank)**2/m)/denominator
            ix, iy = endx, endy
        raw = accum*(total-1)/total
        result[column] = (raw-1)/deviation
    return result


def ad_statistics(x, y, already_sorted=False):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1] or min(len(x),len(y)) < 1:
        raise ValueError('Expected two nonempty cells-by-genes matrices')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Nonfinite AD input')
    if not already_sorted:
        x, y = np.sort(x,axis=0), np.sort(y,axis=0)
    return _ad_sorted_columns(x,y,ad_standard_deviation(len(x),len(y)))


class ADCalibration:
    def __init__(self, path):
        table = np.genfromtxt(Path(path), delimiter=',', names=True)
        self.statistic = np.asarray(table['statistic'])
        p = np.asarray(table['pvalue'])
        if len(p)<10 or not np.isfinite(self.statistic).all() or not np.isfinite(p).all():
            raise ValueError('Incomplete/nonfinite R AD calibration')
        if np.any(np.diff(self.statistic)<=0) or np.any(p<=0) or np.any(p>1) or np.any(np.diff(p)>1e-10):
            raise ValueError('Invalid R AD calibration probabilities/order')
        self.log_p = np.log(p)

    def pvalues(self, statistic):
        values = np.asarray(statistic,dtype=float)
        if np.isnan(values).any() or np.isposinf(values).any():
            raise ValueError('Nonfinite informative AD statistic')
        return np.exp(np.interp(values,self.statistic,self.log_p,left=0.0,right=np.log(1e-300)))


@njit(cache=True)
def _energy_subset_scores(distances, row_sums, total_sum, groups, n, m):
    scores = np.empty(len(groups))
    for i in range(len(groups)):
        selected = groups[i]
        across, within = 0.0, 0.0
        for j in range(n):
            across += row_sums[selected[j]]
            for k in range(n):within += distances[selected[j],selected[k]]
        cross = across-within
        other = total_sum-2*across+within
        scores[i] = max(0.0,2*cross/(n*m)-within/(n*n)-other/(m*m))
    return scores


def energy_permutation_full_control(x,y,permutations=10000,seed=0,max_pairwise_mib=1024):
    """Exact empirical V statistic, all supplied cells, uniform label permutations.

    Selecting the smaller label group is algebraically identical to permuting
    all labels, while avoiding an O(N^2) rescan for every permutation.
    """
    x,y = np.asarray(x,dtype=float),np.asarray(y,dtype=float)
    if x.ndim!=2 or y.ndim!=2 or x.shape[1]!=y.shape[1] or min(len(x),len(y))<1:
        raise ValueError('Invalid energy samples')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or permutations<1:
        raise ValueError('Invalid energy values/permutation count')
    if len(x)>len(y):x,y=y,x
    n,m=len(x),len(y);total=n+m
    if total*total*8>max_pairwise_mib*1024**2:raise MemoryError('Recorded pairwise budget exceeded, no silent subsampling')
    distances=cdist(np.vstack([x,y]),np.vstack([x,y]))
    rows=distances.sum(1);whole=float(rows.sum())
    observed=float(_energy_subset_scores(distances,rows,whole,np.arange(n)[None,:],n,m)[0])
    rng=np.random.default_rng(seed);exceed=0
    for start in range(0,permutations,128):
        groups=np.array([rng.choice(total,n,replace=False) for _ in range(min(128,permutations-start))])
        scores=_energy_subset_scores(distances,rows,whole,groups,n,m)
        exceed+=int(np.sum(scores>=observed-1e-12))
    return observed,(exceed+1)/(permutations+1)


def strong_perturbation_mask(cells, ad_de_genes):
    """Paper inclusion: >=25 cells, >=50 AD-BH DEGs; NO energy-p/q gate.

    The CRISPRi-specific RNA-repression filter is explicitly inapplicable to
    this CRISPRko assay; DNA editing success is not inferred.
    """
    return (np.asarray(cells)>=25)&(np.asarray(ad_de_genes)>=50)


def response_feature_union(ad_statistics_matrix, eligible, raw_mean, normalized_variance, z_defined):
    stats=np.asarray(ad_statistics_matrix);eligible=np.asarray(eligible,dtype=bool)
    mean=np.asarray(raw_mean);variance=np.asarray(normalized_variance);defined=np.asarray(z_defined,dtype=bool)
    if stats.ndim!=2 or eligible.shape!=(len(stats),) or stats.shape[1]!=len(mean):raise ValueError('Feature dimensions disagree')
    expressed=np.flatnonzero((mean>.25)&defined&np.isfinite(variance))
    selected=set()
    if len(expressed):
        cutoff=np.quantile(variance[defined&np.isfinite(variance)],.7)
        selected.update(expressed[variance[expressed]>=cutoff].tolist())
    valid=np.flatnonzero(defined)
    for row in stats[eligible]:
        finite=valid[np.isfinite(row[valid])]
        selected.update(finite[np.argsort(row[finite],kind='stable')[-min(10,len(finite)):]].tolist())
    return np.asarray(sorted(selected),dtype=int)


def paper_hdbscan(profiles,eligible,features):
    from hdbscan import HDBSCAN
    profiles=np.asarray(profiles,dtype=float);indices=np.flatnonzero(eligible)
    labels=np.full(len(profiles),-1,dtype=int)
    if len(indices)<4:return labels,indices,'insufficient_eligible_profiles'
    values=profiles[indices][:,features]
    if values.shape[1]<2:return labels,indices,'undefined_profile_correlation'
    valid=np.ptp(values,axis=1)>0
    excluded=int((~valid).sum());indices=indices[valid];values=values[valid]
    if len(indices)<4:return labels,indices,'insufficient_defined_profiles'
    distances=pdist(values,metric='correlation')
    if not np.isfinite(distances).all():raise ValueError('Nonfinite profile distance')
    from scipy.spatial.distance import squareform
    labels[indices]=HDBSCAN(metric='precomputed',min_cluster_size=4,min_samples=1,cluster_selection_method='eom').fit_predict(squareform(np.clip(distances,0,2)))
    order=indices[leaves_list(linkage(np.clip(distances,0,2),method='average'))]
    return labels,order,'hdbscan_evaluated' if excluded==0 else 'hdbscan_evaluated_undefined_profiles_excluded'
