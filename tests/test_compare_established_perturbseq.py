import importlib.util
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse

spec=importlib.util.spec_from_file_location('comparison',Path(__file__).parents[1]/'scripts'/'compare_established_perturbseq.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_count_identity_includes_content_and_order():
    x=sparse.csr_matrix(np.array([[0,2,1],[3,0,4]],dtype=np.int32))
    assert m.canonical_counts(x)==m.canonical_counts(x.astype(float))
    assert m.canonical_counts(x)!=m.canonical_counts(x[:,[1,0,2]])
    with pytest.raises(ValueError):m.canonical_counts(x.astype(float)*.5)


def test_barcodes_must_match_uniquely():
    assert m.checked_index(['b','a'],['a','b']).tolist()==[1,0]
    with pytest.raises(ValueError):m.checked_index(['a','a'],['a'])
    with pytest.raises(ValueError):m.checked_index(['a'],['c'])


def test_ntc_z_is_control_sample_standard_deviation():
    x=sparse.csr_matrix([[2.,1.,0],[4.,1.,0],[30.,2.,5.]])
    mean,sd,valid=m.ntc_z_parameters(x,np.array([True,True,False]))
    np.testing.assert_allclose(mean,[3,1,0])
    np.testing.assert_allclose(sd,[np.sqrt(2),0,0])
    assert valid.tolist()==[True,False,False]


def test_reference_excludes_whole_ntc_guide():
    x=np.array([[2.,4.,1.],[6.,8.,1.],[20.,30.,1.],[4.,6.,1.]])
    guides=np.array(['N1','N1','N2','T'])
    table,effect=m.profiles(x,[('N2',[2]),('T',[3])],np.array([1,1,1,0],bool),guides)
    np.testing.assert_allclose(effect[0],[16,24,0])
    np.testing.assert_allclose(effect[1],x[3]-x[:3].mean(0))
    assert table.reference_cells.tolist()==[2,3]


def test_mean_log_and_log_mean_are_distinct():
    x=np.array([[0.,0.,0.],[8.,2.,4.],[1.,1.,1.],[1.,1.,1.]])
    guides=np.array(['T','T','N','N'])
    _,a=m.profiles(x,[('T',[0,1])],np.array([0,0,1,1],bool),guides,mode='lfc')
    _,b=m.profiles(np.log1p(x),[('T',[0,1])],np.array([0,0,1,1],bool),guides)
    assert not np.allclose(a*np.log(2),b)


def test_nearest_controls_whole_guide_exclusion():
    e=np.arange(7,dtype=float)[:,None]
    guides=np.array(['N1','N1','N2','N2','N3','N3','T'])
    nc=m.nearest_controls(e,np.array([1,1,1,1,1,1,0],bool),guides,k=2)
    assert not (guides[nc]==guides[:,None]).any()
    assert nc[-1].tolist()==[5,4]


def test_zero_variance_not_fabricated_correlation():
    assert np.isnan(m.pearson(np.zeros(3),np.arange(3)))
    assert m.pearson(np.arange(4),-np.arange(4))==pytest.approx(-1)
