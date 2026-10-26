import importlib.util
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


spec=importlib.util.spec_from_file_location('public_audit',Path(__file__).parents[1]/'scripts/audit_public_perturbseq.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


def write_counts(h,x):
    x=sparse.csc_matrix(x,dtype=np.float32)
    c=h.create_group('counts');c.attrs['encoding-type']='csc_matrix';c.attrs['shape']=x.shape
    for name in ['data','indices','indptr']:c.create_dataset(name,data=getattr(x,name))
    return c


def test_streamed_counts_preserve_every_entry_and_exact_totals(tmp_path):
    expected=np.array([[0,3,1],[2,0,4],[7,2,0],[0,1,9]])
    with h5py.File(tmp_path/'x.h5','w') as h:
        c=write_counts(h,expected)
        data,indices,ptr,totals,detected,digests=audit.read_counts(c,lambda *_args,**_kw:None,chunk_entries=2)
    got=sparse.csc_matrix((data,indices,ptr),shape=expected.shape).toarray()
    np.testing.assert_array_equal(got,expected)
    np.testing.assert_array_equal(totals,expected.sum(1))
    np.testing.assert_array_equal(detected,(expected>0).sum(1))
    assert set(digests)=={'data','indices'}


def test_fractional_counts_cannot_enter_a_raw_count_analysis(tmp_path):
    with h5py.File(tmp_path/'x.h5','w') as h:
        c=write_counts(h,[[1,.25],[0,3]])
        with pytest.raises(ValueError,match='raw integer'):audit.read_counts(c,lambda *_args,**_kw:None)


def test_duplicate_sparse_entries_are_not_double_counted_as_genes(tmp_path):
    with h5py.File(tmp_path/'x.h5','w') as h:
        c=write_counts(h,[[1,2],[3,0]])
        c['indices'][1]=c['indices'][0]
        with pytest.raises(ValueError,match='duplicate-free'):audit.read_counts(c,lambda *_args,**_kw:None)


def metadata_fixture():
    rows=[]
    for batch in ['A','B','C']:
        for i in range(12):
            rows.extend(dict(batch=batch,gRNA=f'NTC{i}',gene_target='NTC') for _ in range(60))
        for i in range(4):
            rows.extend(dict(batch=batch,gRNA=f'T{i}',gene_target=f'G{i//2}') for _ in range(60))
    return pd.DataFrame(rows)


def test_reference_whole_guides_and_both_test_arms_never_overlap():
    obs=metadata_fixture();table,full,reference,plans,a,b=audit.make_plan(obs,[25],17)
    assert set(a).isdisjoint(b)
    assert not np.any((reference>=0)&(plans[25]['assignment']>=0))
    arm_counts=np.bincount(plans[25]['assignment'][plans[25]['assignment']>=0])
    np.testing.assert_array_equal(arm_counts,np.full(len(arm_counts),25))
    assert table.groupby('group').cells.first().sum()==len(obs)
    again=audit.make_plan(obs,[25],17)
    np.testing.assert_array_equal(plans[25]['assignment'],again[3][25]['assignment'])


def test_guide_target_conflicts_are_rejected():
    obs=metadata_fixture();obs.loc[0,'gene_target']='WRONG'
    with pytest.raises(ValueError,match='multiple targets'):audit.make_plan(obs,[25],17)


def test_sparse_group_normalization_matches_direct_dense_arithmetic():
    x=np.array([[0.,3.],[2.,1.],[6.,2.],[1.,0.]])
    totals=x.sum(1);assign=np.array([0,0,1,1]);scale=10
    expected=(x*(scale/totals)[:,None]).reshape(2,2,2).sum(1)
    for col in range(2):
        rr=np.flatnonzero(x[:,col]);v=x[rr,col]*scale/totals[rr]
        np.testing.assert_allclose(audit.aggregate_column(v,rr,assign,2),expected[:,col])


def test_panel_read_is_exact_and_does_not_read_other_columns(tmp_path):
    x=np.array([[1,0,3,0],[0,2,0,1],[4,1,1,5]],dtype=np.float32)
    with h5py.File(tmp_path/'panel.h5','w') as h:
        c=write_counts(h,x)
        values,rows,ptr,_=audit.read_selected_counts(c,np.array([0,3]),lambda *_args,**_kw:None,1024)
        got=sparse.csc_matrix((values,rows,ptr),shape=(3,2)).toarray()
        np.testing.assert_array_equal(got,x[:,[0,3]])
        with pytest.raises(ValueError,match='read budget'):
            audit.read_selected_counts(c,np.array([0,3]),lambda *_args,**_kw:None,1)


def test_panel_selection_uses_rank_and_coverage_not_effect_flags(tmp_path):
    with h5py.File(tmp_path/'selection.h5','w') as h:
        v=h.create_group('var')
        v.create_dataset('_index',data=np.array(['G0','G1','MT-X','R3','R4'],dtype=h5py.string_dtype()))
        v.create_dataset('highly_variable',data=[True,False,False,True,True])
        v.create_dataset('highly_variable_rank',data=[2.,np.nan,np.nan,1.,0.])
        write_counts(h.create_group('layers'),np.ones((4,5)))
        obs=pd.DataFrame(dict(gene_target=['G0','G1','G1','NTC'],effect=[1e5,-1e5,100,0]))
        selected,mask,targets,payload=audit.select_panel(h,obs,panel_size=2,target_size=1)
        np.testing.assert_array_equal(selected,[1,2,3,4])
        np.testing.assert_array_equal(mask,[False,False,True,True])
        assert targets==['G1'] and payload==4*4*8


def test_predeclared_hvg_panel_does_not_inherit_unrelated_full_rna_filter():
    means=np.array([.01,.1,.4,0.,10.])
    panel=np.array([True,True,True,True,False])
    np.testing.assert_array_equal(audit.correlation_feature_mask(means,panel),[True,True,True,False,False])
    np.testing.assert_array_equal(audit.correlation_feature_mask(means),[False,False,False,False,True])


def test_expression_comparison_uses_disjoint_arms_instead_of_self_correlation():
    entries=pd.DataFrame([dict(batch='A',guide='one',target='G',role='targeting',arm=0,arm_id=0),
                          dict(batch='A',guide='one',target='G',role='targeting',arm=1,arm_id=1)])
    means=np.array([[2.,4.,3.,5.],[3.,2.,6.,4.]])
    table=audit.compare_arms(entries,means,{'A':(np.ones(4),np.ones(4))},{'A':(100,100)},'B_mean_log1p_difference',25,1)
    np.testing.assert_allclose(table.expression_r,audit.correlation(means[0],means[1]))
    assert (table.expression_r<.99).all()


def test_known_same_target_signal_is_recovered_against_unrelated_guides():
    rng=np.random.default_rng(21);signal=rng.normal(size=400);baseline=np.full(400,10.)
    rows=[];profiles=[]
    for name,target,effect in [('g1','G',signal),('g2','G',signal),('g3','H',rng.normal(size=400)),('c1','NTC',np.zeros(400)),('c2','NTC',np.zeros(400))]:
        for arm in [0,1]:
            rows.append(dict(batch='A',guide=name,target=target,role='ntc' if target=='NTC' else 'targeting',arm=arm,arm_id=len(profiles)))
            profiles.append(baseline+effect+rng.normal(scale=.15,size=400))
    result=audit.compare_arms(pd.DataFrame(rows),np.array(profiles),{'A':(baseline,baseline)},{'A':(100,100)},'B_mean_log1p_difference',25,7)
    assert result.loc[result.category.eq('same_target_different_guides'),'r'].min()>.95
    assert result.loc[result.category.eq('ntc_pair_background'),'r'].abs().max()<.2


@pytest.mark.skipif(os.name!='posix',reason='Production launcher uses Linux SIGALRM and resource accounting')
@pytest.mark.parametrize('panel_mode',[False,True])
def test_end_to_end_fixture_preserves_source_and_records_selection_scope(tmp_path,panel_mode):
    obs=metadata_fixture();rng=np.random.default_rng(8);x=rng.poisson(2,size=(len(obs),128)).astype(np.float32)
    p=tmp_path/'source.h5ad'
    with h5py.File(p,'w') as h:
        root=h.create_group('layers');write_counts(root,x)
        o=h.create_group('obs')
        for name in ['batch','gRNA','gene_target']:o.create_dataset(name,data=np.array(obs[name],dtype=h5py.string_dtype()))
        for name,values in dict(_index=np.arange(len(obs)).astype(str),channel=obs.batch,n_gRNA=np.ones(len(obs)),
                                n_UMI_counts=x.sum(1),n_genes=(x>0).sum(1),target_knockdown=np.zeros(len(obs)),
                                exceeds_ntc_fdr=np.ones(len(obs),bool),perturbation_cross_batch_significant=np.ones(len(obs),bool)).items():
            if name in ['_index','channel']:values=np.asarray(values,dtype=h5py.string_dtype())
            o.create_dataset(name,data=values)
        v=h.create_group('var');genes=np.array(['G0','G1',*[f'RNA{i}' for i in range(126)]],dtype=h5py.string_dtype())
        v.create_dataset('_index',data=genes);v.create_dataset('gene_ids',data=genes)
        v.create_dataset('highly_variable',data=np.ones(128,bool))
        v.create_dataset('highly_variable_rank',data=np.arange(128,dtype=float))
    before=p.read_bytes();s=p.stat();out=tmp_path/'result'
    config=dict(output=str(out),source=str(p),expected_bytes=s.st_size,expected_mtime_ns=s.st_mtime_ns,
                arm_sizes=[25],seed=4,historical_sha256='fixture',max_seconds=30)
    if panel_mode:config.update(feature_mode='preexisting_hvg_panel',panel_size=110,target_panel_size=2,max_payload_bytes=4*1024**2)
    cp=tmp_path/'config.json';cp.write_text(json.dumps(config));audit.run(cp)
    assert p.read_bytes()==before
    result=json.loads((out/'summary.json').read_text())
    assert result['scope']=='author_selected_strong_perturbation_cohort_only'
    assert not (out/'owner.lock').exists()
    proof=json.loads((out/'input_proof.json').read_text())
    assert proof['author_total_count_mismatches']==(None if panel_mode else 0)
    assert proof['normalization_totals_source']==('author_obs' if panel_mode else 'recomputed_all_counts')
    if panel_mode:
        with np.load(out/'count_panel_checkpoint.npz') as cache:
            cached=sparse.csc_matrix((cache['data'],cache['indices'],cache['indptr']),shape=(len(obs),110)).toarray()
            np.testing.assert_array_equal(cached,x[:,:110])
        receipt=json.loads((out/'count_panel_checkpoint.json').read_text())
        assert receipt['cells']==len(obs) and receipt['gene_columns']==110
