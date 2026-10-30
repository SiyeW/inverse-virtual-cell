from pathlib import Path
import itertools,sys
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import energy_distance
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from guide_rna_energy import energy_v,balanced_values,unbiased_from_balanced_v,monte_carlo_p,test_pair as run_pair,reference_rows
from diagnose_energy_reference import independent_reference

def test_released_energy_definition_matches_independent_scipy_1D():
    X=np.array([0.,1.,4.,5.]);Y=np.array([2.,3.,7.,9.]);D=cdist(np.r_[X,Y][:,None],np.r_[X,Y][:,None])
    v=energy_v(D,np.arange(8)<4)
    assert np.isclose(v,energy_distance(X,Y)**2)

def test_balanced_batch_and_unbiased_formula_on_all_labelings():
    x=np.array([0.,1.,3.,8.])[:,None];D=cdist(x,x);masks=[]
    for take in itertools.combinations(range(4),2):
        m=np.zeros(4);m[list(take)]=1;masks.append(m)
    masks=np.array(masks);v=balanced_values(D,masks)
    u=unbiased_from_balanced_v(v,D);direct=[]
    for m in masks:
        a=np.flatnonzero(m);b=np.flatnonzero(1-m)
        direct.append(2*D[np.ix_(a,b)].mean()-D[np.ix_(a,a)].sum()/2-D[np.ix_(b,b)].sum()/2)
    assert np.allclose(u,direct) and np.isclose(u.mean(),0)
    assert np.array_equal(np.argsort(v),np.argsort(u))

def test_ties_are_not_false_positive_and_p_never_zero():
    p,b=monte_carlo_p(0,np.zeros(99));assert p==1 and b==99
    p,b=monte_carlo_p(100,np.zeros(99));assert p==.01 and b==0

def test_deterministic_permutations_and_detectable_shift():
    rng=np.random.default_rng(45);X=rng.normal(size=(25,5));Y=rng.normal(size=(25,5))+4
    a,na=run_pair(X,Y,19,199);b,nb=run_pair(X,Y,19,199)
    assert a==b and np.array_equal(na,nb) and a['pvalue']<=.01 and a['energy_unbiased']>0

def test_constant_profiles_have_no_difference():
    r,n=run_pair(np.ones((10,4)),np.ones((10,4)),42,99)
    assert r['pvalue']==1 and np.isclose(r['energy_unbiased'],0)

def test_ntc_query_is_fully_excluded_from_its_reference():
    labels=np.array(['A']*10+['N1']*10+['N2']*10+['N3']*10);roles=np.array(['targeting']*10+['non_targeting_control']*30);barcodes=np.array([str(i) for i in range(40)])
    a=reference_rows('N1',labels,roles,barcodes,10,42)
    assert len(a)==10 and not set(a)&set(range(10,20)) and (roles[a]=='non_targeting_control').all()
    assert np.array_equal(a,reference_rows('N1',labels,roles,barcodes,10,42))
    b=independent_reference('N1',labels,roles,barcodes,10,42)
    assert not set(b)&set(range(10,20)) and (roles[b]=='non_targeting_control').all()
    assert np.array_equal(b,independent_reference('N1',labels,roles,barcodes,10,42))
