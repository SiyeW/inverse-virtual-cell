from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from compare_guide_variation import make_plan,Accumulator,process_public,bind_inputs
from public_count_cache import axis_hash


def fixture():
    gids=['A1','A2','B1','B2','C1','C2','N1','N2'];values=[[1,2,3,4]]*2+[[4,3,2,1]]*2+[[2,4,1,3]]*2+[[2,2,3,3]]*2
    x=np.repeat(np.array(values,float),4,axis=0)
    cells=pd.DataFrame(dict(cell_barcode=['c'+str(i) for i in range(32)],guide_id=np.repeat(gids,4),total_counts=10,
                            role=np.repeat(['targeting']*6+['non_targeting_control']*2,4)))
    records=pd.DataFrame(dict(guide=gids,target=['A','A','B','B','C','C','NTC','NTC'],cells=4,
                             source_control_type=['targeting']*6+['non_targeting_control']*2,direct_feature_index=np.nan))
    pairs=pd.DataFrame([dict(target=g,guide_a=g+'1',guide_b=g+'2',correlation_defined=True) for g in ['A','B','C']])
    support=pd.DataFrame([dict(target=g,source_control_type='targeting',pairs=1) for g in ['A','B','C']])
    genes=np.array(['g1','g2','g3','g4'])
    return cells,records,pairs,support,genes,x


def test_disjoint_equal_groups_and_coverage_only_case():
    cells,records,pairs,support,genes,x=fixture();p=make_plan(cells,records,pairs,support,genes,genes,n=2,repeats=2)
    for (guide,size,repeat),group in p['draws'].groupby(['guide','cells','repeat']):
        a,b=[set(map(int,row.analysis_rows.split(';'))) for row in group.sort_values('side').itertuples()]
        assert len(a)==len(b)==size and not a&b
        assert set(cells.iloc[list(a|b)].guide_id)=={guide}
    assert p['plan'].iloc[p['case_index']].target=='A'
    assert np.allclose(np.asarray(p['operator'].sum(1)),1)


def test_distances_and_reference_cancellation(tmp_path):
    cells,records,pairs,support,genes,x=fixture();p=make_plan(cells,records,pairs,support,genes,genes,n=2,repeats=2)
    means=(p['operator']@sparse.csr_matrix(x)).toarray()
    expected=x[::4];ref=np.array([2,2,3,3])
    acc=Accumulator(p,expected,ref,tmp_path/'checkpoint.npz');acc.consume(np.arange(4),means)
    _,s,cases=acc.finish('toy',genes,genes)
    assert (s.loc[s.comparison.isin(['within_guide','same_target','NTC']),'mean_absolute_difference']<1e-12).all()
    assert (s.loc[s.comparison.eq('different_target'),'mean_absolute_difference']>0).all()
    a,b=means[0],means[4]
    assert np.allclose(np.log2((a+1)/(ref+1))-np.log2((b+1)/(ref+1)),np.log2(a+1)-np.log2(b+1))
    assert cases.cells_per_group.eq(2).all()


def test_checkpoint_reuse_and_input_change_refused(tmp_path):
    cells,records,pairs,support,genes,x=fixture();p=make_plan(cells,records,pairs,support,genes,genes,n=2,repeats=2)
    bind_inputs(p,cells,10,'source')
    means=(p['operator']@sparse.csr_matrix(x)).toarray();path=tmp_path/'x.npz'
    a=Accumulator(p,x[::4],np.array([2,2,3,3]),path);a.consume([0,1],means[:,:2])
    b=Accumulator(p,x[::4],np.array([2,2,3,3]),path)
    assert b.done.tolist()==[True,True,False,False]
    with pytest.raises(ValueError,match='replayed'):b.consume([0,1],means[:,:2])
    bind_inputs(p,cells,20,'source')
    with pytest.raises(ValueError,match='contract'):Accumulator(p,x[::4],np.array([2,2,3,3]),path)


def test_direct_target_mask_same_for_all_comparisons():
    cells,records,pairs,support,genes,x=fixture();records.loc[records.target.eq('A'),'direct_feature_index']=0
    p=make_plan(cells,records,pairs,support,genes,genes,n=2,repeats=2)
    assert not p['keep'][0,0] and p['keep'].shape==(3,4)


def test_public_cache_resume_does_not_reread_completed_source(tmp_path,monkeypatch):
    cells,records,pairs,support,genes,x=fixture();base=tmp_path/'base';tables=base/'tables';tables.mkdir(parents=True);cache=base/'raw_count_cache';cache.mkdir()
    obs=cells.set_index('cell_barcode');obs['source_row']=np.arange(len(obs));obs['main_analysis_stage1']=True
    obs['gene_target']=np.repeat(records.target.to_numpy(dtype=str),4)
    obs.to_csv(tables/'all_released_cells_with_additional_QC.csv.gz')
    pd.DataFrame({'gene_name':genes},index=genes).to_csv(tables/'all_RNA_feature_selection.csv.gz')
    records.to_csv(tables/'guide_effect_metadata.csv',index=False);pairs.to_csv(tables/'within_target_pairs.csv',index=False);support.to_csv(tables/'target_support.csv',index=False)
    np.savez_compressed(tables/'rna_effect_profiles.npz',means=x[::4],reference_means=np.tile([2,2,3,3],(8,1)),guides=records.guide.to_numpy(dtype=str),gene_ids=genes)
    entries=[]
    for first in [0,2]:
        matrix=sparse.csc_matrix(x[:,first:first+2]);p=cache/f'b{first}.npz'
        np.savez_compressed(p,data=matrix.data,rows=matrix.indices,ptr=matrix.indptr)
        entries.append(dict(file=p.name,first=first,end=first+2,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    manifest=cache/'manifest.json';manifest.write_text(json.dumps(dict(contract={'rows':axis_hash(obs.index),'genes':axis_hash(genes)},blocks=entries)))
    cfg=dict(baseline=str(base),cache_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),scale=10)
    out=tmp_path/'result';out.mkdir();original=Path.read_bytes
    def interrupted(p):
        if p==cache/'b2.npz':raise TimeoutError('synthetic interruption')
        return original(p)
    monkeypatch.setattr(Path,'read_bytes',interrupted)
    with pytest.raises(TimeoutError):process_public(cfg,genes,out,2,2,1,lambda *a,**k:None)
    def resumed(p):
        if p==cache/'b0.npz':raise AssertionError('Already accepted raw block was reread')
        return original(p)
    monkeypatch.setattr(Path,'read_bytes',resumed)
    _,acc,_,_,_=process_public(cfg,genes,out,2,2,1,lambda *a,**k:None)
    assert acc.done.all() and acc.max_error<1e-8
