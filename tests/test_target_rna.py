import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import pytest
import assess_target_rna as target_module
from assess_target_rna import resolve_targets,summarize_gene,validate_cache,run_dataset
from compare_agreement_overlay import sha


def test_exact_id_not_inferred_host_or_symbol():
    c=pd.DataFrame(dict(guide=['a','b','c'],target_gene_id=['PB.1','missing','NTC'],target_gene_name=['wrong alias','same','NTC'],source_control_type=['targeting','targeting','non_targeting_control']))
    g=pd.DataFrame(dict(gene_id=['PB.1','id2'],gene_name=['same','PB.1_host']))
    m=resolve_targets(c,g)
    assert m.mapping_status.tolist()==['exact_catalogue_id','missing_from_full_reference','not_targeted']
    assert m.rna_column.tolist()==[0,-1,-1]


def test_zeros_total_normalization_and_pooled_batch_reference():
    cells=pd.DataFrame(dict(guide_id=['NT1','NT2','g','g'],role=['non_targeting_control']*2+['targeting']*2,batch=['A','B','A','A']))
    design=pd.DataFrame(dict(guide=['g'],target_gene_name=['gene'],source_control_type=['targeting']))
    # Same target raw counts in two cells with different full RNA totals.
    raw=np.array([2,4,0,2]);totals=np.array([10000,10000,10000,20000])
    t=summarize_gene(raw*10000/totals,cells,design,'id','gene').iloc[0]
    assert t.mean_per_10000==.5 and t.NTC_mean_per_10000==3
    assert t.RNA_remaining_ratio==1/6 and t.detected_fraction==.5
    assert t.same_batch_RNA_remaining_ratio==.25 and t.NTC_cells==2


def test_zero_reference_and_unobserved_guide_are_not_zero_efficiency():
    cells=pd.DataFrame(dict(guide_id=['NT','g'],role=['non_targeting_control','targeting']))
    cells=pd.concat([cells.iloc[[0]],cells],ignore_index=True)
    design=pd.DataFrame(dict(guide=['g','absent'],target_gene_name=['gene']*2,source_control_type=['targeting']*2))
    r=summarize_gene([0,0,0],cells,design,'id','gene')
    assert r.RNA_remaining_ratio.isna().all() and r.loc[0,'mean_per_10000']==0
    assert r.loc[1,'cells']==0 and np.isnan(r.loc[1,'mean_per_10000'])


def test_target_checkpoint_identity_and_corruption_fail_closed(tmp_path):
    p=tmp_path/'counts';p.write_bytes(b'counts');r=tmp_path/'receipt.json'
    r.write_text(json.dumps(dict(contract='fixed',sha256=sha(p))))
    validate_cache(p,r,'fixed')
    with pytest.raises(ValueError):validate_cache(p,r,'changed')
    p.write_bytes(b'broken')
    with pytest.raises(ValueError):validate_cache(p,r,'fixed')


def test_saved_lab_counts_survive_evaluation_failure_without_source_reread(tmp_path,monkeypatch):
    from scipy import sparse
    c=pd.DataFrame(dict(guide=['g','NT'],target_gene_id=['id','NTC'],target_gene_name=['gene','NTC'],source_control_type=['targeting','non_targeting_control']))
    cat=tmp_path/'catalogue.csv';c.to_csv(cat,index=False)
    cells=pd.DataFrame(dict(cell_barcode=['a','b','c'],guide_id=['g','NT','NT'],role=['targeting','non_targeting_control','non_targeting_control'],total_counts=[10,10,20]))
    genes=pd.DataFrame(dict(gene_id=['id','other'],gene_name=['gene','other']))
    def loader(*_):return sparse.csr_matrix([[0,10],[2,8],[4,16]]),cells,genes,dict(exact=True)
    monkeypatch.setattr(target_module,'load_counts',loader)
    original=target_module.summarize_gene
    monkeypatch.setattr(target_module,'summarize_gene',lambda *_:(_ for _ in ()).throw(ValueError('deliberate post-cache interruption')))
    d=dict(id='fixture',kind='counts',count_input={},target_catalogue=str(cat),target_catalogue_sha256=sha(cat),metadata_sha256={})
    out=tmp_path/'output'
    with pytest.raises(ValueError,match='interruption'):run_dataset(d,out,lambda *_a,**_k:None)
    assert (out/'target_counts.json').exists()
    monkeypatch.setattr(target_module,'load_counts',lambda *_:(_ for _ in ()).throw(AssertionError('Source must not reopen')))
    monkeypatch.setattr(target_module,'summarize_gene',original)
    r=run_dataset(d,out,lambda *_a,**_k:None)
    table=pd.read_csv(out/'guide_target_RNA.csv')
    assert r['exact_id_targets']==1 and table.iloc[0].RNA_remaining_ratio==0


def test_public_cache_exact_id_zero_in_mean_and_completed_block_resume(tmp_path):
    from public_count_cache import axis_hash
    b=tmp_path/'baseline';(b/'tables').mkdir(parents=True);(b/'raw_count_cache').mkdir()
    obs=pd.DataFrame(dict(guide_id=['g_1','g_1','NT_1','NT_1'],gene_target=['gene','gene','NTC','NTC'],batch=['A']*4,source_row=[0,1,2,3],main_analysis_stage1=[True]*4,total_counts=[10,20,10,20]),index=['a','b','c','d'])
    genes=pd.DataFrame(dict(gene_name=['gene','other']),index=pd.Index(['id','other_id'],name='gene_id'))
    obs.to_csv(b/'tables/all_released_cells_with_additional_QC.csv.gz');genes.to_csv(b/'tables/all_RNA_feature_selection.csv.gz')
    block=b/'raw_count_cache/b.npz'
    np.savez(block,data=np.array([2,4,10,18,6,16]),rows=np.array([2,3,0,1,2,3]),ptr=np.array([0,2,6]))
    manifest=dict(contract=dict(rows=axis_hash(obs.index),genes=axis_hash(genes.index)),blocks=[dict(first=0,end=2,file='b.npz',bytes=block.stat().st_size,sha256=sha(block))])
    m=b/'raw_count_cache/manifest.json';m.write_text(json.dumps(manifest))
    c=pd.DataFrame(dict(id=['g-1','NT-1'],target_gene_id=['id','Non-Targeting'],target_gene_name=['gene','Non-Targeting']))
    cat=tmp_path/'catalogue.csv';c.to_csv(cat,index=False)
    d=dict(id='fixture',kind='public_cache',baseline=str(b),cache_manifest_sha256=sha(m),target_catalogue=str(cat),target_catalogue_sha256=sha(cat),metadata_sha256={})
    out=tmp_path/'out';r=run_dataset(d,out,lambda *_a,**_k:None)
    t=pd.read_csv(out/'guide_target_RNA.csv').iloc[0]
    assert t.cells==2 and t.mean_per_10000==0 and t.NTC_mean_per_10000==2000 and t.RNA_remaining_ratio==0
    block.unlink()  # Completed result blocks must not reopen an old source block.
    run_dataset(d,out,lambda *_a,**_k:None)
