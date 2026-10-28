from pathlib import Path
import sys,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pytest
from scipy import stats
from target_rna_precision import guide_interval,matched_ntc_range,saved_target_block,check_original_point,dataset,contracts


def test_interval_matches_standard_scipy_BCa_and_original_mean():
    v=np.array([0,1,0,2,1,3,0,2,1,1,2,4.],float)
    r=guide_interval(v,2,123,resamples=1000)
    expected=stats.bootstrap((v,),np.mean,n_resamples=1000,batch=128,vectorized=True,method='BCa',rng=np.random.default_rng(123))
    assert r['RNA_remaining_ratio']==v.mean()/2
    assert np.allclose([r['guide_ci_low'],r['guide_ci_high']],np.array(expected.confidence_interval)/2)
    assert r['guide_ci_reference_fixed'] and r['guide_ci_status']=='estimated'


def test_no_zero_width_RNA_interval_from_all_zero_or_constant_observations():
    r=guide_interval(np.zeros(10),1,31)
    assert r['RNA_remaining_ratio']==0 and np.isnan(r['guide_ci_high'])
    assert r['guide_ci_status']=='all_zero_no_RNA_mean_interval'
    assert np.isclose(r['zero_detection_fraction_upper95'],stats.beta.ppf(.95,1,10))
    assert guide_interval(np.ones(10),1,31)['guide_ci_status']=='constant_observations_no_interval'


def test_equal_n_NTC_groups_exclude_their_own_cells_from_reference():
    # All six unordered samples of two distinct observations have ratios
    # 3/7,4/6,5/5,5/5,6/4,7/3 against their complementary two cells.
    r=matched_ntc_range([1,2,3,4],2,912,resamples=1000)
    assert r['NTC_null_sample_size']==2 and r['NTC_null_reference_excludes_selected_cells']
    assert np.isclose(r['NTC_equal_n_low'],3/7) and np.isclose(r['NTC_equal_n_high'],7/3)


def test_degenerate_low_expression_NTC_reference_is_not_silently_dropped():
    r=matched_ntc_range([0,0,0,3],2,21,resamples=1000)
    assert r['NTC_equal_n_low']==0 and np.isinf(r['NTC_equal_n_high'])
    assert r['NTC_null_zero_reference_draws']>0
    with pytest.raises(ValueError):matched_ntc_range([0,0,0],1,1)


def test_selected_exact_count_checkpoint_resumes_without_original_payload(tmp_path):
    p=tmp_path/'raw.npz';np.savez(p,data=np.array([1,2,3,4]),rows=np.array([0,2,1,3]),ptr=np.array([0,2,4]))
    entry=dict(first=10,end=12,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    remap=np.array([0,-1,1,2]);dest=tmp_path/'selected.npz'
    m,read=saved_target_block(p,entry,[10,11],3,remap,dest,'fixed')
    assert read and np.array_equal(m.toarray(),[[1,0],[2,0],[0,4]])
    p.unlink()
    again,read=saved_target_block(p,entry,[10,11],3,remap,dest,'fixed')
    assert not read and np.array_equal(m.toarray(),again.toarray())
    with pytest.raises(ValueError):saved_target_block(p,entry,[10,11],3,remap,dest,'changed')


def test_changed_point_or_guide_membership_stops_the_branch():
    from types import SimpleNamespace
    row=SimpleNamespace(cells=3,detected_cells=2,mean_per_10000=1.,NTC_mean_per_10000=2.)
    check_original_point([0,1,2],2,row)
    with pytest.raises(ValueError):check_original_point([0,2,3],2,row)


def test_completed_method_contract_cannot_hide_changed_resample_count():
    a={'seed':1,'resamples':100,'expected_commit':'fixed'}
    b=dict(a,resamples=200)
    assert contracts({'input':'same'},a)[0]==contracts({'input':'same'},b)[0]
    assert contracts({'input':'same'},a)[1]!=contracts({'input':'same'},b)[1]


def test_exact_lab_panel_end_to_end_no_raw_matrix_reader(tmp_path):
    import json
    import pandas as pd
    from scipy import sparse
    from compare_agreement_overlay import sha
    from assess_target_rna import summarize_gene
    base=tmp_path/'baseline';base.mkdir()
    cells=pd.DataFrame(dict(cell_barcode=[f'c{i}' for i in range(30)],guide_id=['g']*10+['NT']*20,role=['targeting']*10+['non_targeting_control']*20,total_counts=[100]*30))
    counts=np.array([0,0,1,0,2,0,1,0,0,1]+[1,2]*10)
    panel=sparse.csc_matrix(counts[:,None]);sparse.save_npz(base/'target_counts.npz',panel)
    cells.to_csv(base/'cells.csv.gz',index=False);pd.DataFrame(dict(gene_id=['id'],gene_name=['gene'])).to_csv(base/'genes.csv.gz',index=False)
    metadata={n:sha(base/n) for n in ['cells.csv.gz','genes.csv.gz']}
    receipt=dict(positions=[0],sha256=sha(base/'target_counts.npz'),metadata=metadata)
    (base/'target_counts.json').write_text(json.dumps(receipt))
    design=pd.DataFrame(dict(guide=['g'],target_gene_name=['gene'],source_control_type=['targeting']))
    table=summarize_gene(counts*100,cells,design,'id','gene');table.insert(0,'dataset','fixture');table.to_csv(base/'guide_target_RNA.csv',index=False)
    d=dict(id='fixture',kind='counts',target_result=str(base),point_table_sha256=sha(base/'guide_target_RNA.csv'),target_counts_sha256=sha(base/'target_counts.npz'),target_counts_receipt_sha256=sha(base/'target_counts.json'),cache_metadata_sha256=metadata,expected_cells=30,expected_NTC_cells=20)
    cfg=dict(seed=1,resamples=100,expected_commit='fixture')
    out=tmp_path/'out';result=dataset(d,cfg,out,lambda *_a,**_k:None)
    actual=pd.read_csv(out/'guide_target_RNA_precision.csv')
    assert result['original_points_unchanged'] and actual.iloc[0].guide_ci_status=='estimated'
    assert np.isclose(actual.iloc[0].RNA_remaining_ratio,table.iloc[0].RNA_remaining_ratio)
    (base/'genes.csv.gz').write_bytes(b'changed')
    with pytest.raises(ValueError,match='metadata changed'):dataset(d,cfg,out,lambda *_a,**_k:None)
