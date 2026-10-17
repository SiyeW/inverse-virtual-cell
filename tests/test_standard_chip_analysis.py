from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from perturbation_statistics import energy_permutation, energy_permutation_batched
from run_chip01_standard import main_mask, normalize_counts, energy_representation


def test_batched_energy_matches_original_same_permutations():
    rng=np.random.default_rng(12);x=rng.normal(size=(12,3));y=rng.normal(size=(15,3))
    old=energy_permutation(x,y,permutations=199,seed=4)
    new=energy_permutation_batched(x,y,permutations=199,seed=4,batch_size=17)
    np.testing.assert_allclose(old,new,rtol=1e-12)
    assert energy_permutation_batched(np.zeros((12,2)),np.zeros((15,2)),99)==(0,1)
    assert energy_permutation_batched(np.zeros((12,2)),np.ones((15,2))*10,99)[1]==.01


def test_policy_is_explicit_and_does_not_mutate():
    obs=pd.DataFrame(dict(qc_pass_stage1=[True,True,False,True],guide_status=['single_guide','multi_guide','single_guide','single_guide'],main_analysis_stage1=[True,False,False,True]))
    original=obs.copy(deep=True)
    assert main_mask(obs).tolist()==[True,False,False,True]
    pd.testing.assert_frame_equal(obs,original)
    obs.loc[1,'main_analysis_stage1']=True
    with pytest.raises(ValueError):main_mask(obs)


def test_normalization_and_energy_exclude_custom_features_and_preserve_counts():
    rng=np.random.default_rng(5);counts=sparse.csr_matrix(rng.poisson(2,size=(60,40)))
    original=counts.copy();totals=np.asarray(counts.sum(1)).ravel()
    normalized=normalize_counts(counts,totals)
    np.testing.assert_allclose(np.asarray(normalized.sum(1)).ravel(),10000,rtol=1e-6)
    var=pd.DataFrame(dict(gene_id=['ENSG'+str(i) for i in range(38)]+['PB.x','PB.y'],gene_name=['g'+str(i) for i in range(40)]))
    scores,indices,variance,_=energy_representation(counts,totals,var,np.arange(60)<30)
    assert scores.shape==(60,20) and max(indices)<38 and np.isfinite(scores).all()
    assert len(variance)==20 and (counts!=original).nnz==0
