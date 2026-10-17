from pathlib import Path
import sys,warnings
import numpy as np
import pytest
from scipy.stats import anderson_ksamp
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from replogle2022_statistics import (ad_statistics,ad_standard_deviation,energy_permutation_full_control,
    _energy_subset_scores,strong_perturbation_mask,response_feature_union,paper_hdbscan)


@pytest.mark.parametrize('sizes',[(1,30),(2,2),(3,7),(17,31),(20,400)])
@pytest.mark.parametrize('tied',[False,True])
def test_ad_matches_scipy_midrank(sizes,tied):
    rng=np.random.default_rng(14);x=rng.normal(size=(sizes[0],5));y=rng.normal(size=(sizes[1],5))
    if tied:x=np.round(x);y=np.round(y)
    actual=ad_statistics(x,y)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        expected=[anderson_ksamp([x[:,g],y[:,g]]).statistic for g in range(5)]
    np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=1e-12)
    np.testing.assert_allclose(ad_statistics(y,x),actual,rtol=2e-12,atol=1e-12)


def test_ad_constant_and_invalid():
    assert np.isneginf(ad_statistics(np.zeros((10,1)),np.zeros((20,1)))[0])
    with pytest.raises(ValueError):ad_standard_deviation(1,1)
    with pytest.raises(ValueError):ad_statistics(np.array([[np.nan],[1]]),np.ones((3,1)))


def test_energy_subset_is_exact_v_statistic():
    rng=np.random.default_rng(11);a=rng.normal(size=(18,4));d=cdist(a,a)
    groups=np.array([rng.choice(18,5,False) for _ in range(10)])
    result=_energy_subset_scores(d,d.sum(1),d.sum(),groups,5,13)
    for selected,got in zip(groups,result):
        other=np.setdiff1d(np.arange(18),selected)
        expected=2*d[np.ix_(selected,other)].mean()-d[np.ix_(selected,selected)].mean()-d[np.ix_(other,other)].mean()
        np.testing.assert_allclose(got,expected,atol=2e-14)


def test_energy_null_and_full_control_not_silently_capped():
    x=np.zeros((20,2));y=np.zeros((301,2))
    assert energy_permutation_full_control(x,y,99)==(0,1)
    effect,p=energy_permutation_full_control(x,y+10,99)
    assert effect>20 and p==.01
    assert energy_permutation_full_control(x,y+10,99)==energy_permutation_full_control(x,y+10,99)
    with pytest.raises(MemoryError):energy_permutation_full_control(x,y,max_pairwise_mib=.001)


def test_paper_admission_is_not_energy_gated_and_feature_union():
    np.testing.assert_array_equal(strong_perturbation_mask([25,24,50,80],[50,100,49,60]),[True,False,False,True])
    stats=np.tile(np.arange(40),(4,1));mean=np.ones(40);var=np.arange(40,dtype=float)
    features=response_feature_union(stats,[True,False,False,False],mean,var,np.ones(40,bool))
    assert set(range(30,40)).issubset(features)
    assert len(features)>=12


def test_paper_hdbscan_and_no_fit_when_empty():
    rng=np.random.default_rng(4)
    profile=np.vstack([np.tile([1,1,0,0],(5,1)),np.tile([0,0,1,1],(5,1))])+rng.normal(scale=.001,size=(10,4))
    labels,order,status=paper_hdbscan(profile,np.ones(10,bool),np.arange(4))
    assert status=='hdbscan_evaluated' and len(set(labels))==2 and -1 not in labels
    labels,order,status=paper_hdbscan(profile,np.zeros(10,bool),np.arange(4))
    assert status=='insufficient_eligible_profiles' and len(order)==0 and np.all(labels==-1)


def test_core_csv_identity_ignores_storage_dtype_but_rejects_real_changes():
    import pandas as pd
    from run_chip01_literature import validate_core_policy_identity
    original=pd.DataFrame(dict(source_row=[1,3],total_counts=[10,20],is_control=[True,False],
                          guide_id=['NTC1','g1'],perturbed_gene=['ctrl','GENE1']),
                          index=pd.Index(pd.array(['AA-1','BB-1'],dtype='string')))
    policy=original.copy();policy.index=pd.Index(np.array(['AA-1','BB-1'],dtype=object),dtype=object)
    validate_core_policy_identity(policy,original)
    with pytest.raises(ValueError):validate_core_policy_identity(policy.iloc[::-1],original)
    broken=policy.copy();broken.loc['BB-1','guide_id']='other'
    with pytest.raises(ValueError):validate_core_policy_identity(broken,original)
    broken=policy.copy();broken.loc['BB-1','total_counts']=21
    with pytest.raises(ValueError):validate_core_policy_identity(broken,original)
