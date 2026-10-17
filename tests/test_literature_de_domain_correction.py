from pathlib import Path
import argparse,json,sys
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from correct_literature_de_domain import run
from analyze_expression_structure import sha256


def test_full_domain_correction_keeps_source_and_energy_but_changes_bh_family(tmp_path):
    source=tmp_path/'parent';source.mkdir();destination=tmp_path/'corrected'
    domain=np.array([False,True,True,True]);names=np.array(['A','B','C','D'])
    pd.DataFrame(dict(gene_name=names,relative_z_defined=domain),index=names).to_csv(source/'gene_statistics.csv.gz')
    index=[]
    for i,(target,cells) in enumerate([('P1',12),('P2',6)]):
        tested=cells>=10;p=np.array([1e-217,.01,.04,.2])
        pd.DataFrame(dict(gene_id=names,gene_name=names,response_z_defined=domain,
                         ad_statistic=[500,5,3,1],ad_test_informative=True,ad_pvalue=p,ad_bh_q=p,
                         mw_tested=tested,mw_pvalue=p if tested else np.ones(4),mw_bh_q=p,
                         mean_normalized=[1,2,3,4])).to_csv(source/f'target_{i:03d}_de.csv.gz',index=False)
        index.append(dict(target=target,file=f'target_{i:03d}_de.csv.gz',cells=cells,tested_genes=4,
                          ad_de_genes=3,mw_tested=tested,mw_de_genes=3 if tested else np.nan))
    pd.DataFrame(index).to_csv(source/'de_summary.csv',index=False)
    np.savez_compressed(source/'response_profiles.npz',targets=np.array(['P1','P2']),values=np.ones((2,4)),
                        ad_statistics=np.ones((2,4)),z_defined=domain,raw_mean=np.ones(4),raw_variance=np.ones(4))
    for name in ['cell_policy.csv.gz','energy_representation.npz','energy_results.csv','energy_progress.csv','direct_target_masking.csv']:(source/name).write_bytes(b'unchanged verified artifact')
    (source/'frozen_parameters.json').write_text('{}');(source/'summary.json').write_text(json.dumps(dict(rna_features=4,scientific_result=False)))
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in source.iterdir()}
    (source/'output_manifest.json').write_text(json.dumps(manifest));(source/'status.json').write_text(json.dumps(dict(stage='complete',code_commit='fixture-parent')))
    original_hashes={p.name:sha256(p) for p in source.iterdir()}
    run(argparse.Namespace(source_effects=source,output=destination))
    table=pd.read_csv(destination/'target_000_de.csv.gz')
    assert table.de_gene_testable.tolist()==domain.tolist()
    assert np.isnan(table.loc[0,['ad_pvalue','ad_bh_q','mw_pvalue','mw_bh_q']].to_numpy(dtype=float)).all()
    np.testing.assert_allclose(table.loc[domain,'ad_bh_q'],[.03,.06,.2])
    second=pd.read_csv(destination/'target_001_de.csv.gz');assert not second.mw_tested.any() and second.mw_pvalue.isna().all()
    assert pd.read_csv(destination/'de_summary.csv').tested_genes.tolist()==[3,3]
    assert json.loads((destination/'summary.json').read_text())['tested_rna_features']==3
    assert np.isneginf(np.load(destination/'response_profiles.npz')['ad_statistics'][:,0]).all()
    assert {p.name:sha256(p) for p in source.iterdir()}==original_hashes
    assert sha256(destination/'energy_results.csv')==sha256(source/'energy_results.csv')
