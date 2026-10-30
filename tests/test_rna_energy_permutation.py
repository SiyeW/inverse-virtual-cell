import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from rna_energy_permutation import subset_statistics,permutation_samples,depth_strata,permutation_pvalue


@pytest.mark.parametrize('nc,nq,n',[(17,0,3),(11,4,4),(9,6,5)])
def test_unequal_size_and_partition_blocks(nc,nq,n):
    rng=np.random.default_rng(4)
    pc=rng.normal(size=(nc+nq,5));d=cdist(pc,pc);rows=d.sum(axis=1);total=float(rows.sum())
    masks=permutation_samples(len(pc),np.arange(n),31,902)
    got=subset_statistics(masks,d[:nc,:nc],d[nc:,:nc],d[nc:,nc:],rows,total)
    for index,selected in enumerate(masks):
        other=np.setdiff1d(np.arange(len(pc)),selected);a=d[np.ix_(selected,selected)].sum();b=d[np.ix_(other,other)].sum();cross=d[np.ix_(selected,other)].mean()
        assert got[index,0]==pytest.approx(2*cross-a/n**2-b/len(other)**2,abs=1e-11)
        assert got[index,1]==pytest.approx(2*cross-a/(n*(n-1))-b/(len(other)*(len(other)-1)),abs=1e-11)


def test_restricted_permutations_preserve_covariate_counts_without_dropping_cells():
    values=np.linspace(1000,30000,60);blocks=depth_strata(values,5);observed=np.array([0,1,2,15,27,40,59])
    draws=permutation_samples(60,observed,71,93,blocks)
    assert np.array_equal(draws,permutation_samples(60,observed,71,93,blocks))
    for draw in draws:
        assert len(np.unique(draw))==len(observed)
        assert np.array_equal(np.bincount(blocks[draw],minlength=5),np.bincount(blocks[observed],minlength=5))


def test_ties_cannot_become_significant():
    assert permutation_pvalue(0,np.zeros(499))==(1.,499)
    assert permutation_pvalue(10,np.zeros(499))==(1/500,0)


def test_invalid_groups_and_covariates_fail():
    with pytest.raises(ValueError):permutation_samples(10,np.array([1,1]),10,0)
    with pytest.raises(ValueError):depth_strata([1,0,2])
