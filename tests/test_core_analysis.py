"""New-chip portability, historical numerical parity and integrity regression tests."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))

import perturbseq_core_effects as core
from perturbseq_core_io import load_config, read_catalogue, read_source, select_features, sha256
from perturbseq_core_report import verify_output
from run_core_analysis import run
from test_chip01_report_alignment import fixture


def inputs(tmp_path, chip='chip07', missing_optional_roles=False):
    counts, obs, overlay, _ = fixture()
    if missing_optional_roles:
        keep = ~obs.guide_id.isin(['ig','p1','p2']).to_numpy()
        counts, obs = counts[keep], obs.iloc[np.flatnonzero(keep)].copy()
        overlay = overlay.loc[~overlay.guide.isin(['ig','p1','p2'])].copy()
    obs.total_counts = np.asarray(counts.sum(1)).ravel()
    obs['chip_id'] = chip
    obs['biological_replicate'] = 'not_established'
    var = pd.DataFrame({'gene_name':list('PTABCD'),'gene_id':[f'feature{i}' for i in range(6)]}, index=[f'feature{i}' for i in range(6)])
    source = tmp_path/f'{chip}.h5ad'
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(source)
    catalogue = tmp_path/f'{chip}_catalogue.csv'
    overlay[['guide','target_gene_name','source_control_type']].to_csv(catalogue,index=False)
    config = dict(schema_version=1, run_id=chip+'_fixture', chip_id=chip,
        expected_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source=dict(path=str(source),bytes=source.stat().st_size,shape=list(counts.shape),sha256=sha256(source)),
        catalogue=dict(path=str(catalogue),sha256=sha256(catalogue)),
        control_target='Non-Targeting', features=dict(mode='pooled_raw_mean',raw_mean_gt=.5),
        normalization_target_umi='median_ntc',min_guide_cells=10,matched_cells=10,matched_repeats=2,
        seed=20261022,threads=1,smoke_prefix_rows=1000,allowed_output_parent=str(tmp_path),output=str(tmp_path/'run'))
    path = tmp_path/'config.json'
    path.write_text(json.dumps(config),encoding='utf8')
    return path, config, counts, obs, overlay, var


def test_extracted_effects_are_numerically_identical_to_accepted_implementation():
    from chip01_report_alignment import effect_profiles, target_support
    counts, obs, overlay, _ = fixture()
    reference = obs.guide_id.isin(['n1','n2']).to_numpy()
    args = (counts,obs,overlay,reference,list('PTABCD'),list('PTABCD'),10,10,3,20261022)
    old, new = effect_profiles(*args), core.effect_profiles(*args)
    pd.testing.assert_frame_equal(old[0],new[0])
    for i in range(1,5):np.testing.assert_array_equal(old[i],new[i])
    pd.testing.assert_frame_equal(old[5],new[5]);assert old[6]==new[6]
    a,b=target_support(counts,obs,overlay,old[0],old[3],reference)
    c,d=core.target_support(counts,obs,overlay,new[0],new[3],reference)
    pd.testing.assert_frame_equal(a,c);pd.testing.assert_frame_equal(b,d)


def test_reader_preserves_original_counts_rows_metadata_and_does_one_formal_hash(tmp_path):
    _,config,counts,obs,_,var=inputs(tmp_path)
    saved=Path(config['source']['path']).read_bytes()
    x,selected,all_obs,features,proof=read_source(config['source'])
    np.testing.assert_array_equal(x.toarray(),counts.toarray())
    assert selected.index.tolist()==obs.index.tolist() and proof['full_hash_verified']
    assert features.index.tolist()==var.index.tolist()
    assert Path(config['source']['path']).read_bytes()==saved
    bad=dict(config['source'],sha256='0'*64)
    with pytest.raises(ValueError,match='SHA256'):read_source(bad)


@pytest.mark.parametrize('problem',['duplicate','missing','wrong_target','wrong_role'])
def test_catalogue_conflicts_stop_without_guessing(tmp_path,problem):
    _,config,_,obs,_,_=inputs(tmp_path)
    p=Path(config['catalogue']['path']);table=pd.read_csv(p)
    if problem=='duplicate':table=pd.concat([table,table.iloc[:1]])
    elif problem=='missing':table=table.loc[table.guide.ne('n1')]
    elif problem=='wrong_target':table.loc[table.guide.eq('t1'),'target_gene_name']='WRONG'
    else:table.loc[table.guide.eq('n1'),'source_control_type']='targeting'
    table.to_csv(p,index=False)
    with pytest.raises(ValueError):read_catalogue(dict(path=str(p),sha256=sha256(p)),obs,'Non-Targeting')


def test_frozen_ids_preserve_exact_order_and_reject_ambiguity(tmp_path):
    _,_,counts,_,_,var=inputs(tmp_path)
    path=tmp_path/'features.csv'
    pd.DataFrame(dict(gene_id=var.index[[0,2,5]])).to_csv(path,index=False)
    spec=dict(mode='frozen_ids',path=str(path),sha256=sha256(path))
    np.testing.assert_array_equal(select_features(counts,var,spec),[0,2,5])
    pd.DataFrame(dict(gene_id=var.index[[5,0]])).to_csv(path,index=False)
    spec['sha256']=sha256(path)
    with pytest.raises(ValueError,match='order'):select_features(counts,var,spec)


@pytest.mark.parametrize('missing_roles',[False,True])
def test_new_chip_end_to_end_report_atlas_and_read_only_acceptance(tmp_path,missing_roles):
    path,config,counts,_,_,_=inputs(tmp_path,missing_optional_roles=missing_roles)
    source_before=sha256(config['source']['path'])
    summary=run(path,smoke=True)
    assert summary['chip_id']=='chip07' and not summary['scientific_result']
    assert summary['analysis_cells']==counts.shape[0] and summary['inferential_de']=='not_run'
    result=verify_output(Path(config['output']))
    assert result['accepted'] and not result['raw_source_reopened']
    text=(Path(config['output'])/'分析报告.md').read_text(encoding='utf8')
    assert 'chip07' in text and 'Day 0/Day 10' in text and '工程前缀测试' in text
    assert 'figures/05_matched_rna_effects.png' in text
    assert sha256(config['source']['path'])==source_before
    assert not (Path(config['output'])/'null_runtime').exists()
    with pytest.raises(FileExistsError):run(path,smoke=True)
    bad=Path(config['output'])/'figures/05_matched_rna_effects.png'
    with bad.open('ab') as stream:stream.write(b'corrupt')
    with pytest.raises(ValueError,match='Result differs'):verify_output(Path(config['output']))


def test_zero_pair_targets_are_not_silently_dropped():
    x,obs,overlay,_=fixture()
    reference=obs.guide_id.isin(['n1','n2']).to_numpy()
    records,_,_,effects,_,_,_=core.effect_profiles(x,obs,overlay,reference,list('PTABCD'),list('PTABCD'),14,10,2,0)
    targets,pairs=core.target_support(x,obs,overlay,records,effects,reference)
    assert len(targets)==3 and targets.eligible_guides.le(1).all() and pairs.empty
    assert 'pearson' in pairs and targets.assessment.eq('insufficient_two_guide_profiles').all()


def test_overall_only_diagnostic_retains_scope_tables_and_valid_report(tmp_path):
    path,config,counts,_,_,_=inputs(tmp_path,chip='company_merged')
    config.update(max_target_cards=0, input_scope='company_merged_records_diagnostic',
                  report_notice='公司合并矩阵的诊断；未确认六次捕获的细胞来源。')
    path.write_text(json.dumps(config),encoding='utf8')
    summary=run(path,smoke=True)
    output=Path(config['output'])
    assert summary['target_cards']==0 and summary['max_target_cards']==0
    assert summary['input_scope']=='company_merged_records_diagnostic'
    assert not list((output/'figures/targets').glob('*.png'))
    assert len(pd.read_csv(output/'tables/target_support.csv'))==3
    assert pd.read_csv(output/'tables/target_figure_index.csv').empty
    assert config['report_notice'] in (output/'分析报告.md').read_text(encoding='utf8')
    assert '未生成逐靶点图片' in (output/'分析报告.md').read_text(encoding='utf8')
    assert verify_output(output)['accepted']


@pytest.mark.parametrize('cards',[-1,True,1.5,'all'])
def test_invalid_card_limit_is_rejected(tmp_path,cards):
    path,config,*_=inputs(tmp_path)
    config['max_target_cards']=cards
    path.write_text(json.dumps(config),encoding='utf8')
    with pytest.raises(ValueError,match='max_target_cards'):load_config(path)


def test_merged_diagnostic_requires_visible_notice(tmp_path):
    path,config,*_=inputs(tmp_path)
    config['input_scope']='company_merged_records_diagnostic'
    path.write_text(json.dumps(config),encoding='utf8')
    with pytest.raises(ValueError,match='report_notice'):load_config(path)


@pytest.mark.parametrize('notice',[None,'',' ',False,15])
def test_merged_diagnostic_rejects_empty_or_nontext_notice(tmp_path,notice):
    path,config,*_=inputs(tmp_path)
    config.update(input_scope='company_merged_records_diagnostic',report_notice=notice)
    path.write_text(json.dumps(config),encoding='utf8')
    with pytest.raises(ValueError,match='report_notice'):load_config(path)


def test_standalone_cli_has_no_legacy_import_dependency():
    code="import sys;sys.path.insert(0,'scripts');import run_core_analysis;assert not any(x.startswith(('run_chip01','chip01_','replogle2022','perturbation_statistics')) for x in sys.modules)"
    subprocess.run([sys.executable,'-c',code],cwd=ROOT,check=True,capture_output=True,text=True)
    result=subprocess.run([sys.executable,'scripts/run_core_analysis.py','--help'],cwd=ROOT,check=True,capture_output=True,text=True)
    assert 'verify' in result.stdout


def test_archive_manifest_is_byte_exact_and_all_old_paths_are_retired():
    manifest=json.loads((ROOT/'archive/20261023_legacy_chip_analysis/archive_manifest.json').read_text())
    assert len(manifest['files'])==19
    for entry in manifest['files']:
        path=ROOT/entry['archived_path']
        assert path.stat().st_size==entry['bytes'] and sha256(path)==entry['sha256']
        assert not (ROOT/entry['former_path']).exists()


def test_unsafe_config_and_string_booleans_rejected(tmp_path):
    path,config,_,obs,_,var=inputs(tmp_path)
    config['chip_id']='../chip01'
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError,match='identifier'):load_config(path)
    from perturbseq_core_io import validate_metadata
    obs['is_control']=obs.is_control.astype(str)
    with pytest.raises(ValueError,match='boolean'):validate_metadata(obs,var)


def test_preparation_names_follow_chip_and_protect_raw_sources(tmp_path,monkeypatch):
    import run_core_analysis as runner
    source=tmp_path/'matrix.h5'
    matrix=sparse.csc_matrix(np.vstack([np.full((200,4),4),np.ones((1,4)),[[7,0,9,2]]]).astype('int32'))
    with h5py.File(source,'w') as handle:
        group=handle.create_group('matrix')
        for name,value in [('data',matrix.data),('indices',matrix.indices),('indptr',matrix.indptr),('shape',matrix.shape)]:group[name]=value
        group['barcodes']=np.asarray(['a','b','c','d'],dtype='S')
        features=group.create_group('features')
        features['id']=np.asarray([f'ENSG{i}' for i in range(201)]+['GUIDE'],dtype='S')
        features['name']=np.asarray([f'GENE{i}' for i in range(200)]+['MT-CO1','guide1'],dtype='S')
        features['feature_type']=np.asarray(['Gene Expression']*201+['CRISPR Guide Capture'],dtype='S')
    calls=tmp_path/'calls.csv'
    pd.DataFrame(dict(cell_barcode=['a','c','d'],num_features=[1,2,1],
        feature_call=['NTC','NTC|PB-guide','PB-guide'],num_umis=['7','8|1','2'])).to_csv(calls,index=False)
    ref=tmp_path/'ref.csv'
    pd.DataFrame(dict(name=['NTC','PB-guide'],target_gene_name=['Non-Targeting','PB.custom'])).to_csv(ref,index=False)
    metrics=tmp_path/'metrics.csv';pd.DataFrame(dict(metric=[1])).to_csv(metrics,index=False)
    config=dict(chip_id='chip09',expected_commit='test',matrix_h5=str(source),protospacer_calls=str(calls),
        feature_reference=str(ref),metrics_summary=str(metrics),input_hashes_json='fixture-only',
        input_sha256={k:sha256(p) for k,p in zip(['matrix_h5','protospacer_calls','feature_reference','metrics_summary'],[source,calls,ref,metrics])},
        allowed_output_parent=str(tmp_path),output_root=str(tmp_path/'baseline'))
    path=tmp_path/'prepare.json';path.write_text(json.dumps(config))
    monkeypatch.setattr(runner,'revision',lambda *args,**kwargs:dict(commit='fixture',working_tree_dirty=False))
    before=sha256(source)
    runner.prepare(path)
    output=Path(config['output_root'])
    assert (output/'data/chip09_all_cells_stage1_qc.h5ad').is_file()
    assert not (output/'data/chip01_all_cells_stage1_qc.h5ad').exists()
    assert sha256(source)==before and verify_output(output)['accepted']
    with pytest.raises(FileExistsError):runner.prepare(path)
