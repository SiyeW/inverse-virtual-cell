from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from chip01_report_alignment import (TYPES,annotation_overlay,reference_masks,
    leave_guide_out_reference,expression_effect,effect_profiles,pair_effects,
    target_support,qc_flow,unique_direct_feature)


def fixture():
    guide=np.repeat(['n1','n2','ig','p1','p2','t1','t2'],[12,13,10,12,11,13,14])
    target={'n1':'ctrl','n2':'ctrl','ig':'ctrl','p1':'P','p2':'P','t1':'T','t2':'T','zero':'Q'}
    kind={'n1':TYPES[0],'n2':TYPES[0],'ig':TYPES[1],'p1':TYPES[2],'p2':TYPES[2],'t1':TYPES[3],'t2':TYPES[3],'zero':TYPES[3]}
    obs=pd.DataFrame(dict(guide_id=guide,perturbed_gene=[target[g] for g in guide],
        is_control=np.isin(guide,['n1','n2','ig']),source_row=np.arange(len(guide)),
        guide_status='single_guide',qc_pass_stage1=True,main_analysis_stage1=True,
        total_counts=100.,n_genes_by_counts=50,pct_counts_mt=3.,raw_guide_umis=7,guide_umis=6),index=[f'c{i}' for i in range(len(guide))])
    lib=[];annotations=[]
    for g in sorted(target):
        coarse='NTC' if kind[g] in TYPES[:2] else 'targeting'
        name='Non-Targeting' if coarse=='NTC' else target[g]
        n=int((guide==g).sum())
        lib.append(dict(name=g,id=g,target_gene_name=name,guide_type=coarse,single_guide_before_qc=n,single_guide_after_stage1_qc=n))
        annotations.append(dict(guide=g,current_target=name,current_type=coarse,comparator_type=kind[g]))
    overlay=annotation_overlay(pd.DataFrame(lib),annotations)
    rng=np.random.default_rng(8);matrix=sparse.csr_matrix(rng.poisson(3,(len(obs),6)).astype(float))
    core=np.isin(guide,['n1','n2','ig'])
    return matrix,obs,overlay,core


def test_complete_overlay_retains_zero_and_original_coarse_roles():
    _,_,overlay,_=fixture()
    assert len(overlay)==8 and overlay.loc[overlay.guide.eq('zero'),'single_guide_after_stage1_qc'].iloc[0]==0
    assert overlay.loc[overlay.guide.eq('ig'),'guide_type'].iloc[0]=='NTC'
    assert not overlay.source_role_validated_in_this_assay.any()
    rows=overlay.rename(columns={'guide':'name'})
    ann=[dict(guide=r.name,current_target=r.target_gene_name,current_type=r.guide_type,comparator_type=r.source_control_type) for r in rows.itertuples()]
    with pytest.raises(ValueError):annotation_overlay(rows,ann[:-1])
    bad=[dict(x) for x in ann];bad[0]['current_target']='wrong'
    with pytest.raises(ValueError):annotation_overlay(rows,bad)


def test_pure_reference_excludes_intergenic_and_retains_original_mask():
    _,obs,overlay,old=fixture();before=old.copy()
    types,pure,new=reference_masks(obs,overlay,old)
    assert old.sum()==35 and new.sum()==25 and pure.sum()==25
    assert not new[types==TYPES[1]].any()
    np.testing.assert_array_equal(old,before)


def test_reference_is_leave_whole_guide_out_without_source_mutation():
    x,obs,overlay,old=fixture();saved=x.copy()
    _,_,ref=reference_masks(obs,overlay,old)
    mean,indices=leave_guide_out_reference(x,obs.guide_id,ref,'n1')
    assert len(indices)==13 and not obs.iloc[indices].guide_id.eq('n1').any()
    np.testing.assert_allclose(mean,x[obs.guide_id.eq('n2').to_numpy()].toarray().mean(0))
    assert (x!=saved).nnz==0


def test_expression_effect_has_declared_pseudocount_and_signed_zero():
    np.testing.assert_allclose(expression_effect([0,3,1],[0,1,3]),[0,1,-1])
    with pytest.raises(ValueError):expression_effect([0],[0],0)
    with pytest.raises(ValueError):expression_effect([np.nan],[0])
    assert unique_direct_feature('G',['G','H'],['G','G','H']) is None


def test_matched_profiles_use_equal_n_and_nonoverlapping_pure_reference():
    x,obs,overlay,old=fixture();_,_,ref=reference_masks(obs,overlay,old)
    names=['P','T','A','B','C','D']
    result=effect_profiles(x,obs,overlay,ref,names,names,10,10,3,53)
    records,means,centers,effects,first,metrics,splits=result
    assert len(records)==7 and metrics.cells.eq(10).all() and metrics.reference_overlap_cells.eq(0).all()
    assert records.loc[records.guide.eq('n1'),'reference_cells'].iloc[0]==13
    assert records.loc[records.guide.eq('p1'),'reference_cells'].iloc[0]==25
    assert np.isnan(effects[records.guide.eq('p1'),0]).all()
    for row in splits:
        assert len(row['source_rows'])==len(set(row['source_rows']))==10
        selected=obs.iloc[row['source_rows']]
        assert selected.guide_id.eq(row['guide']).all()
    again=effect_profiles(x,obs,overlay,ref,names,names,10,10,3,53)
    np.testing.assert_array_equal(first,again[4]);assert splits==again[6]
    for i in range(len(records)):
        actual=expression_effect(means[i],centers[i]);keep=np.isfinite(effects[i])
        np.testing.assert_allclose(actual[keep],effects[i,keep])


def test_target_support_preserves_unassessed_and_direction_is_descriptive():
    x,obs,overlay,old=fixture();_,_,ref=reference_masks(obs,overlay,old)
    records,_,_,effects,_,_,_=effect_profiles(x,obs,overlay,ref,list('ABCDEF'),list('ABCDEF'),10,10,2,0)
    targets,pairs=target_support(x,obs,overlay,records,effects,ref)
    assert len(targets)==3 and len(pairs)==2
    assert targets.loc[targets.target.eq('Q'),'assessment'].iloc[0]=='insufficient_two_guide_profiles'
    assert targets.loc[targets.target.eq('P'),'cells_without_largest'].iloc[0]==11
    assert 'pvalue' not in pairs and 'hit' not in targets
    assert not pair_effects(np.zeros(3),np.ones(3))['correlation_defined']


def test_cell_flow_conserves_input_and_never_reselects():
    _,obs,_,_=fixture();before=obs.copy(deep=True)
    obs.loc['c0','guide_status']='multi_guide';obs.loc['c0','main_analysis_stage1']=False
    flow=qc_flow(obs)
    assert flow.total_cells.sum()==len(obs)
    assert flow.loc[flow.guide_status.eq('multi_guide'),'total_cells'].iloc[0]==1
    bad=before.copy();bad.loc['c0','main_analysis_stage1']=False
    with pytest.raises(ValueError):qc_flow(bad)


def test_bounded_h5ad_smoke_reads_contiguous_counts_and_metadata(tmp_path):
    import anndata as ad
    from run_chip01_report_alignment import read_prefix_smoke
    x,obs,_,_=fixture()
    obs['total_counts']=np.asarray(x.sum(1)).ravel()
    var=pd.DataFrame(index=list('ABCDEF'))
    source=tmp_path/'source.h5ad';ad.AnnData(X=x,obs=obs,var=var).write_h5ad(source)
    before=source.read_bytes()
    actual,meta,identity=read_prefix_smoke(source,obs,var,20)
    np.testing.assert_array_equal(actual.toarray(),x[:20].toarray())
    assert meta.index.tolist()==obs.index[:20].tolist()
    assert identity['source_prefix_rows']==20 and not identity['source_full_hash_recomputed']
    assert source.read_bytes()==before


def test_static_output_smoke_keeps_fixed_cards_and_relative_links(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    from plot_chip01_report_alignment import technical_figures,reference_figure,effect_background_figures,target_atlas
    from run_chip01_report_alignment import example_means,role_summary
    x,obs,overlay,old=fixture();_,_,ref=reference_masks(obs,overlay,old)
    records,means,centers,effects,first,metrics,_=effect_profiles(x,obs,overlay,ref,list('ABCDEF'),list('ABCDEF'),10,10,2,0)
    for folder in ['figures/targets','tables']:(tmp_path/folder).mkdir(parents=True,exist_ok=True)
    technical_figures(obs,overlay,qc_flow(obs),tmp_path)
    impact=pd.DataFrame(dict(historical_mean=[1,2,3],pure_core_mean=[1.1,2.1,3.1],reference_log2_change=[.02,.03,.04],all_ntc_vs_pure_core_log2_change=[.04,.02,.03]))
    reference_figure(impact,35,25,tmp_path)
    examples,selections=example_means(x,obs,overlay,ref,10,8)
    assert all(s['reference_overlap_cells']==0 for s in selections)
    effect_background_figures(records,effects,first,metrics,examples,tmp_path)
    targets,pairs=target_support(x,obs,overlay,records,effects,ref)
    cards=target_atlas(overlay,records,effects,targets,pairs,tmp_path)
    assert len(cards)==2 and all((tmp_path/r['figure']).exists() for r in cards)
    assert (tmp_path/'逐靶点支持图册.md').is_file()
    assert len(list((tmp_path/'figures').glob('*.png')))==6
    roles=role_summary(overlay,records,tmp_path)
    assert sum(r['library_guides'] for r in roles)==8
