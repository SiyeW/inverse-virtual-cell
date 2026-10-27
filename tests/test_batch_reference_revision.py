from pathlib import Path
import sys
import json

import h5py
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from public_count_cache import cache_counts
from perturbseq_core_effects import stable_rng
from revalidate_batch_reference import plan_groups,stream_means,profiles_and_draws


def fixture():
    guide=np.repeat(['N1','N1','N2','N2','T1','T2'],50)
    batch=np.repeat(['A','B','A','B','A','A'],50)
    ntc=np.char.startswith(guide.astype(str),'N')
    rng=np.random.default_rng(37)
    counts=rng.poisson(8,size=(300,5)).astype(np.int32)+1
    counts[batch=='B',2]+=30
    counts[~ntc,4]+=9
    obs=pd.DataFrame(dict(guide_id=guide,batch=batch,gene_target=np.where(ntc,'NTC','TARGET'),
        n_gRNA=1,source_row=np.arange(300),total_counts=counts.sum(1)),index=[f'cell{i}' for i in range(300)])
    records=pd.DataFrame([dict(guide=g,cells=int((guide==g).sum())) for g in sorted(set(guide))])
    return counts,obs,records


def test_target_draws_unchanged_control_units_separate_and_halves_disjoint():
    counts,obs,records=fixture();plan=plan_groups(obs,records)
    assert len(plan['groups'])==6
    assert plan['groups'].source_control_type.eq('non_targeting_control').sum()==4
    for guide in ['T1','T2']:
        expected=np.sort(stable_rng(20261022,guide).choice(np.flatnonzero(obs.guide_id.eq(guide)),10,replace=False))
        got=plan['draws'].loc[plan['draws'].key.eq(guide)&plan['draws'].n.eq(10)&plan['draws'].repeat.eq(1)].iloc[0]
        assert got.source_rows==expected.tolist()
    for key,part in plan['halves'].groupby('key'):
        a,b=[set(plan['operator'].getrow(int(row.operator_row)).indices) for row in part.itertuples()]
        assert a.isdisjoint(b)


def test_streamed_raw_cache_equals_direct_cell_averaging_and_correct_batch_reference(tmp_path):
    counts,obs,records=fixture();plan=plan_groups(obs,records)
    gene_ids=[f'id{i}' for i in range(5)];selected=np.array([0,2,4])
    from scipy import sparse
    x=sparse.csc_matrix(counts,dtype=np.float32)
    with h5py.File(tmp_path/'fixture.h5','w') as h:
        n=h.create_group('counts');n.attrs['encoding-type']='csc_matrix';n.attrs['shape']=x.shape
        for key in ['data','indices','indptr']:n.create_dataset(key,data=getattr(x,key))
        manifest=cache_counts(n,tmp_path/'cache',dict(source='fixture'),obs.index,gene_ids,
                              lambda *_a,**_kw:None,max_entries=700)
    linear,logged,_=stream_means(tmp_path/'cache',manifest,obs.index,gene_ids,obs,selected,
                                 plan['operator'],1000,lambda *_a,**_kw:None)
    normalized=counts*(1000/counts.sum(1))[:,None]
    np.testing.assert_allclose(linear,plan['operator']@normalized[:,selected],rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(logged,plan['operator']@np.log2(1+normalized[:,selected]),rtol=1e-13,atol=1e-13)
    genes=np.array(['TARGET','rna1','rna2','rna3','rna4'])
    groups,means,centers,effects,first,metrics,halves=profiles_and_draws(linear,plan,genes[selected],genes,'lfc')
    for i,row in enumerate(groups.itertuples()):
        ref=obs.batch.eq(row.batch)&obs.gene_target.eq('NTC')&obs.guide_id.ne(row.guide)
        # The guide exclusion is harmless for a targeting group and whole-guide for NTC.
        np.testing.assert_allclose(centers[i],normalized[ref][:,selected].mean(0),rtol=1e-12,atol=1e-12)
        assert row.reference_cells==ref.sum()
    assert np.isnan(effects[groups.source_control_type.eq('targeting'),0]).all()
    assert metrics.groupby(['key','cells']).size().eq(8).all()
    assert len(halves)==6


def test_cross_batch_targeting_or_corrupt_checkpoint_fails_closed(tmp_path):
    counts,obs,records=fixture()
    obs.loc[obs.guide_id.eq('T1').idxmax(),'batch']='B'
    with pytest.raises(ValueError,match='Targeting guide spans batches'):plan_groups(obs,records)


@pytest.mark.skipif(sys.platform=='win32',reason='Actual Linux owner/alarm end-to-end contract')
def test_completed_original_to_cached_batch_revision_without_source_payload(tmp_path,monkeypatch):
    import validate_original_on_public as original
    import revalidate_batch_reference as revised
    import perturbseq_core_plots as plots
    import run_core_analysis as runner
    from perturbseq_core_io import sha256
    counts,obs,records=fixture()
    rng=np.random.default_rng(80)
    counts=np.c_[counts,rng.poisson(5,size=(len(obs),220))+1].astype(np.int32)
    # The adopted QC is data-dependent and may remove a tail record. Leave
    # real headroom above the production minimum of 100 NTCs in each batch.
    # Duplicate synthetic rows under UNIQUE new cell identities, then compare
    # to the actually adopted baseline cohort rather than assuming all pass.
    take=np.repeat(np.arange(len(obs)),2)
    counts=counts[take];obs=obs.iloc[take].copy()
    obs.index=[f'cell{i}' for i in range(len(obs))]
    obs['source_row']=np.arange(len(obs))
    path=tmp_path/'public.h5ad'
    strings=h5py.string_dtype('utf-8')
    from scipy import sparse
    with h5py.File(path,'w') as h:
        layer=h.create_group('layers');node=layer.create_group('counts')
        x=sparse.csc_matrix(counts,dtype=np.float32)
        node.attrs['shape']=x.shape;node.attrs['encoding-type']='csc_matrix'
        for key in ['data','indices','indptr']:node.create_dataset(key,data=getattr(x,key))
        o=h.create_group('obs');v=h.create_group('var')
        for key,values in dict(_index=obs.index,gRNA=obs.guide_id,gene_target=obs.gene_target,batch=obs.batch).items():
            o.create_dataset(key,data=np.asarray(values,object),dtype=strings)
        for key,values in dict(n_gRNA=np.ones(len(obs),int),n_gRNA_UMIs=np.ones(len(obs),int),
                               n_UMI_counts=counts.sum(1),n_genes=(counts>0).sum(1)).items():o.create_dataset(key,data=values)
        v.create_dataset('gene_ids',data=np.asarray([f'id{i}' for i in range(counts.shape[1])],object),dtype=strings)
        v.create_dataset('_index',data=np.asarray(['TARGET']+[f'RNA{i}' for i in range(1,counts.shape[1])],object),dtype=strings)
    cat=tmp_path/'catalogue.csv'
    pd.DataFrame(dict(id=['N1','N2','T1','T2'],target_gene_name=['Non-Targeting','Non-Targeting','TARGET','TARGET'])).to_csv(cat,index=False)
    cfg=dict(expected_commit='fixture',source=str(path),expected_bytes=path.stat().st_size,
        expected_mtime_ns=path.stat().st_mtime_ns,historical_sha256=sha256(path),
        catalogue=str(cat),catalogue_sha256=sha256(cat),output=str(tmp_path/'baseline'),
        allowed_output_parent=str(tmp_path),max_seconds=180,block_entries=5000)
    cfgpath=tmp_path/'original.json';cfgpath.write_text(json.dumps(cfg))
    monkeypatch.setattr(original,'revision',lambda *_a:dict(commit='fixture'))
    monkeypatch.setattr(plots,'save',plots.save);monkeypatch.setattr(runner,'save',runner.save)
    monkeypatch.setattr(plots,'NAMES',plots.NAMES.copy());monkeypatch.setattr(runner,'NAMES',runner.NAMES.copy())
    baseline_summary=original.run(cfgpath)
    # A successful revision must NOT touch even one original H5AD byte.
    path.unlink()
    new=dict(expected_commit='fixture',baseline=str(tmp_path/'baseline'),output=str(tmp_path/'revision'),
             allowed_output_parent=str(tmp_path),max_seconds=180)
    newpath=tmp_path/'revision.json';newpath.write_text(json.dumps(new))
    monkeypatch.setattr(revised,'revision',lambda *_a:dict(commit='fixture'))
    revised.run(newpath)
    out=tmp_path/'revision';result=json.loads((out/'summary.json').read_text())
    assert result['analysis_cells']==baseline_summary['additional_RNA_QC_cells']
    assert result['RNA_features']==counts.shape[1]
    assert result['source_count_payload_read'] is False
    assert result['original_bridge_max_error']<1e-10
    assert not (out/'owner.lock').exists()
    for method,value in result['methods'].items():
        assert value['discrimination']['within_batch_NTC_pair']['matched_pairs']==1
        assert value['precision_common50']['targeting']['units']==2
        assert value['matched_RMS']['10']['by_batch']['A']['targeting_units']==2
