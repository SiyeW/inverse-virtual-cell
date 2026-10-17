from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.stats import hypergeom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from finish_chip01_analysis import read_gmt,overrepresentation,profile_clusters


def test_gmt_supports_official_and_enrichr_formats(tmp_path):
    p=tmp_path/'a.gmt'
    p.write_text('Official\tR-HSA-1\tA\tB\nEnrichr R-HSA-2\t\tB\tC\nMouse\tR-MMU-3\tA\n')
    sets=read_gmt(p)
    assert set(sets)=={'R-HSA-1','R-HSA-2'} and sets['R-HSA-2'][1]=={'B','C'}


def test_ora_uses_tested_background_and_all_pathways_for_bh():
    sets={'R-HSA-1':('one',{'A','B','outside'}),'R-HSA-2':('two',{'C','D'})}
    table=overrepresentation({'A','B','outside'},set('ABCDE'),sets,min_size=1)
    first=table.set_index('pathway_id').loc['R-HSA-1']
    assert first.query_genes==2 and first.overlap==2 and first.background_genes==5
    assert first.pvalue==hypergeom.sf(1,5,2,2) and first.bh_q==min(1,2*first.pvalue)
    empty=overrepresentation(set(),set('ABCDE'),sets,min_size=1)
    assert (empty.pvalue==1).all() and (empty.bh_q==1).all()


def test_profile_groups_are_similarity_not_forced_all_members():
    rng=np.random.default_rng(1);signal=np.arange(50)
    values=np.vstack([signal+rng.normal(0,.1,50) for _ in range(3)]+[rng.normal(size=50),np.zeros(50)])
    labels,order,cols=profile_clusters(values,np.ones(5,dtype=bool))
    assert labels[0]==labels[1]==labels[2]>=0 and labels[3]==labels[4]==-1
    assert len(order)==4 and len(cols)>1
    labels,_,_=profile_clusters(values,np.zeros(5,dtype=bool))
    assert (labels==-1).all()
    labels,_,cols=profile_clusters(np.tile(signal,(3,1)),np.ones(3,dtype=bool))
    assert (labels==-1).all() and len(cols)==0


def test_functional_workflow_writes_auditable_report_and_preserves_inputs(tmp_path):
    import argparse,json,hashlib
    import pandas as pd
    from finish_chip01_analysis import run
    effects=tmp_path/'effects';structure=tmp_path/'structure';effects.mkdir();structure.mkdir()
    genes=['ENSG'+str(i) for i in range(30)];symbols=['G'+str(i) for i in range(30)]
    targets=['testA','testB','testC']
    universe=pd.DataFrame({'gene_name':symbols},index=pd.Index(genes,name='gene_id'))
    universe.to_csv(effects/'tested_gene_universe.csv')
    rows=[];profiles=[]
    rng=np.random.default_rng(4)
    for i,target in enumerate(targets):
        filename=f'target_{i:03d}_de.csv.gz'
        fc=np.r_[np.ones(13),np.zeros(17)]+rng.normal(0,.001,30)
        pd.DataFrame({'bh_q':np.r_[np.ones(13)*.001,np.ones(17)],'log2_fold_change':fc,
          'ntc_sensitivity_bh_q':np.r_[np.ones(13)*.001,np.ones(17)],'ntc_sensitivity_log2_fold_change':fc},index=genes).to_csv(effects/filename)
        rows.append(dict(target=target,file=filename,cells=60,exploratory_de_genes=13))
        profiles.append(fc)
    pd.DataFrame(rows).to_csv(effects/'de_summary.csv',index=False)
    pd.DataFrame({'target':targets,'cells':[60]*3,'energy_full':[1.,2.,3.],'bh_q':[.001]*3,'ntc_sensitivity_bh_q':[.001]*3}).to_csv(effects/'energy_results.csv',index=False)
    np.savez_compressed(effects/'perturbation_expression_profiles.npz',targets=np.array(targets),values=np.array(profiles))
    (structure/'qc_policy.json').write_text(json.dumps(dict(main_cells=120,controls=30)))
    for folder in [effects,structure]:
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in folder.iterdir()}
        (folder/'output_manifest.json').write_text(json.dumps(manifest))
        (folder/'status.json').write_text(json.dumps(dict(stage='complete',code_commit='fixture')))
    gmt=tmp_path/'toy.gmt'
    gmt.write_text('First\tR-HSA-1\t'+'\t'.join(symbols[:15])+'\nSecond\tR-HSA-2\t'+'\t'.join(symbols[15:])+'\n')
    (tmp_path/'source.json').write_text(json.dumps(dict(sha256=hashlib.sha256(gmt.read_bytes()).hexdigest())))
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in effects.iterdir()}
    output=tmp_path/'functional'
    run(argparse.Namespace(effects=effects,structure=structure,gmt=gmt,output=output))
    assert json.loads((output/'status.json').read_text())['stage']=='complete'
    assert (output/'04_standard_analysis_report.md').exists()
    assert len(pd.read_csv(output/'pathway_enrichment_hits.csv.gz'))>0
    assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in effects.iterdir()}
    from review_chip01_results import run as review
    # The reviewer must describe positive fixture results rather than hard-code
    # the negative outcome observed in the real first-chip analysis.
    (tmp_path/'structure_v2').mkdir()
    import shutil
    for p in structure.iterdir():shutil.copy2(p,tmp_path/'structure_v2'/p.name)
    (tmp_path/'functional_analysis').mkdir()
    for p in output.iterdir():shutil.copy2(p,tmp_path/'functional_analysis'/p.name)
    review(tmp_path,tmp_path/'review')
    assert (tmp_path/'review/descriptive_response_similarity.png').exists()
    assert '入口为空' not in (tmp_path/'review/05_interpretation.md').read_text(encoding='utf8')
