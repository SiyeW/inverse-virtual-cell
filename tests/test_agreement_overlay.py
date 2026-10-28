from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from compare_agreement_overlay import linear_and_log_means,centers_from_means,apply_mask,evaluate,members


def fixture():
    cells=pd.DataFrame(dict(guide_id=['A','A','N','N','tiny','tiny'],role=['targeting','targeting','non_targeting_control','non_targeting_control','non_targeting_control','non_targeting_control'],total_counts=[10]*6))
    rec=pd.DataFrame(dict(guide=['A','N'],target=['T','NTC'],cells=[2,2],source_control_type=['targeting','non_targeting_control'],direct_feature_index=[np.nan,np.nan]))
    return cells,rec


def test_mean_order_and_all_ntc_including_ineligible_guide():
    cells,rec=fixture();x=sparse.csr_matrix([[0,10],[10,0],[2,8],[4,6],[9,1],[9,1]])
    lm,lc,bm,bc,n=linear_and_log_means(x,cells,rec,10)
    assert n==4 and np.array_equal(lm,[[5,5],[3,7]])
    assert np.allclose(lc,[6,4])
    refs=centers_from_means(lm,rec,lc,n)
    assert np.allclose(refs,[[6,4],[9,1]])
    assert np.allclose(bm[0],[np.log2(11)/2,np.log2(11)/2])
    assert not np.allclose(bm[0],np.log2(lm[0]+1))
    assert np.array_equal(x.toarray(),[[0,10],[10,0],[2,8],[4,6],[9,1],[9,1]])


def test_normalization_uses_full_RNA_total_not_panel_sum():
    cells,rec=fixture();x=sparse.csr_matrix([[1],[2],[3],[4],[5],[6]])
    lm,_,_,_,_=linear_and_log_means(x,cells,rec,20)
    assert np.allclose(lm[:,0],[3,7])


def test_changed_cell_identity_fails():
    cells,rec=fixture();rec.loc[0,'cells']=3
    with pytest.raises(ValueError,match='cell count'):linear_and_log_means(sparse.csr_matrix(np.ones((6,2))),cells,rec,10)


def test_direct_target_mask_and_frozen_target_coverage():
    rec=pd.DataFrame(dict(guide=['A','B','C','D','N1','N2'],direct_feature_index=[0,0,np.nan,np.nan,np.nan,np.nan]))
    values=np.array([[50,1,3,2],[2,1,3,2],[1,2,3,4],[4,3,2,1],[2,1,0,3],[0,1,2,3]],float)
    masked=apply_mask(values,rec);assert np.isnan(masked[:2,0]).all() and values[0,0]==50
    plan=pd.DataFrame([dict(target='T',guide_a='A',guide_b='B',other_a='C',other_b='D',NTC_a='N1',NTC_b='N2')])
    support=pd.DataFrame([dict(target='T',source_control_type='targeting',pairs=1,largest_guide_fraction=.6)])
    points,pairs=evaluate(rec,masked,plan,support,'new')
    assert points.median_pair_r.iloc[0]==pytest.approx(1)
    assert set(pairs.comparison)=={'same_target','different_target','NTC'}
    support.loc[0,'target']='wrong'
    with pytest.raises(ValueError,match='coverage'):evaluate(rec,masked,plan,support,'new')


def test_selective_npz_member_is_exact_and_nonobject(tmp_path):
    path=tmp_path/'x.npz';np.savez_compressed(path,means=np.arange(12).reshape(3,4),unused=np.ones((100,100)))
    data,receipt=members(path,['means'])
    assert list(data)==['means'] and np.array_equal(data['means'],np.arange(12).reshape(3,4))
    assert len(receipt['consumed_members'])==1 and receipt['consumed_members'][0]['crc_verified']
