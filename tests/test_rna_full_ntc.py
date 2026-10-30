from pathlib import Path
import sys,itertools
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pytest
from sklearn.metrics import pairwise_distances
from explore_rna_full_ntc import FullControl,energy_u,permutation_variance,moments


def direct(D,idx):
    mask=np.ones(len(D),bool);mask[list(idx)]=False;a=np.asarray(idx);b=np.flatnonzero(mask)
    return energy_u(len(a),len(b),D[np.ix_(a,a)].sum(),D[np.ix_(a,b)].sum(),D[np.ix_(b,b)].sum())


def test_exact_random_label_variance_matches_all_unequal_allocations():
    rng=np.random.default_rng(152);D=pairwise_distances(rng.normal(size=(9,4)))
    T,S2,R2=moments(D)
    for n in [2,3,4,6,7]:
        values=np.array([direct(D,a) for a in itertools.combinations(range(9),n)])
        assert abs(values.mean())<1e-12
        assert np.isclose(permutation_variance(9,n,T,S2,R2),values.var(),rtol=1e-10,atol=1e-12)


def test_full_control_streaming_matches_direct_for_target_and_own_NTC():
    rng=np.random.default_rng(25);nc=rng.normal(size=(17,3));q=rng.normal(size=(5,3));control=FullControl(nc,working_memory=.001)
    D=pairwise_distances(np.vstack([nc,q]));v=control.evaluate(q)
    assert v['reference_records']==17 and v['query_records']==5
    assert np.isclose(v['energy_unbiased'],direct(D,np.arange(17,22)),atol=1e-12)
    T,S2,R2=moments(D);assert np.isclose(v['random_label_SD']**2,permutation_variance(22,5,T,S2,R2),rtol=1e-9)
    own=np.array([1,4,9,13]);v=control.evaluate(nc[own],own)
    assert v['reference_records']==13 and np.isclose(v['energy_unbiased'],direct(pairwise_distances(nc),own),atol=1e-12)
    with pytest.raises(ValueError):control.evaluate(nc[own[:3]],own)


def test_constant_data_has_zero_effect_and_zero_random_variance():
    control=FullControl(np.ones((12,3)),working_memory=.001);r=control.evaluate(np.ones((4,3)))
    assert r['energy_unbiased']==0 and r['random_label_SD']==0


def test_large_pool_equal_off_diagonal_has_zero_variance():
    N=19342
    for n in [10,50,200,532]:
        assert abs(permutation_variance(N,n,N*(N-1),N*(N-1),N*(N-1)**2))<1e-8
