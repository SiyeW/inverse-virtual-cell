"""Complete the frozen RNA-state comparison as descriptive diagnostics.

The NTC review guard is retained and reported, not changed or passed. All
statistic, samples, seeds and permutations match the original method. No
reliable perturbation-response candidate list is produced.
"""
from pathlib import Path
import argparse,hashlib,json,os,socket,sys,time,atexit
import numpy as np
import pandas as pd
from guide_rna_energy import sha,digest,atomic,reference_rows,test_pair


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--extension-commit',required=True);args=parser.parse_args()
    cfg=json.loads(args.config.read_text());task=json.loads(Path(cfg['task_manifest']).read_text());op=Path(task['native_directory']);root=Path(task['canonical_output_root'])
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    assert (op/'run.exit').read_text().strip()=='0' and not (op/'active.lock').exists() and not (op/'reference_diagnosis.lock').exists()
    diagnostic=json.loads((root/'NTC_diagnostic.json').read_text());assert diagnostic['review_required']
    if not (root/'reference_sensitivity/summary.json').exists():raise ValueError('Required NTC sensitivity not completed')
    if sha(Path(__file__).with_name('guide_rna_energy.py'))!=cfg['script_sha256']:raise ValueError('Frozen primary scientific source differs')
    for p,h in {**cfg['input_SHA256'],**cfg['module_SHA256']}.items():
        if sha(p)!=h:raise ValueError('Scientific input or implementation differs')
    lock=op/'diagnostic_extension.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    s=Path('/proc/self/stat').read_text();ticks=int(s[s.rfind(')')+2:].split()[19]);os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,extension_commit=args.extension_commit)).encode());os.close(fd)
    atexit.register(lambda:lock.unlink(missing_ok=True))
    def progress(stage,**kw):atomic(op/'diagnostic_extension.status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kw))
    try:
        with np.load(cfg['input_path'],allow_pickle=False) as z:pc=z['PCA'];bc=z['cell_barcode'];labels=z['guide'];roles=z['role'];totals=z['total_RNA'];RNA_ids=z['RNA_ids']
        assert pc.shape==(71664,40) and len(set(bc))==71664
        assert hashlib.sha256(json.dumps(RNA_ids.tolist()).encode()).hexdigest()==task['intended_RNA_inventory_SHA256']
        cells=pd.read_csv(cfg['all_cells_path'],dtype={'cell_barcode':str,'guide_id':str});counts=cells.groupby('guide_id').size();eligible=cells[cells.guide_id.map(counts).ge(10)]
        assert eligible.cell_barcode.tolist()==bc.tolist() and eligible.guide_id.tolist()==labels.tolist()
        catalogue=pd.read_csv(cfg['catalogue_path'],keep_default_na=False).set_index('guide',verify_integrity=True)
        groups=[]
        for guide in np.unique(labels):
            rows=np.flatnonzero(labels==guide);rr=np.unique(roles[rows]);assert len(rr)==1 and rr[0]==catalogue.loc[guide,'source_control_type']
            groups.append(dict(guide=str(guide),role=str(rr[0]),rows=rows,n=len(rows)))
        groups.sort(key=lambda g:(g['role']!='non_targeting_control',hashlib.sha256(g['guide'].encode()).hexdigest()));assert len(groups)==785
        import pertpy as pt
        from statsmodels.stats.multitest import multipletests
        official=pt.tl.Distance('edistance',obsm_key='X_pca',cell_wise_metric='euclidean').metric_fct
        contract=digest({k:v for k,v in cfg.items() if k!='script_sha256'});results=[]
        for i,g in enumerate(groups):
            file=root/'guides'/f'{i:04d}.json';null=file.with_suffix('.npz');ref=reference_rows(g['guide'],labels,roles,bc,g['n'],cfg['reference_seed'])
            sample_sha=digest(dict(query=bc[g['rows']].tolist(),reference=bc[ref].tolist()))
            if file.exists():
                row=json.loads(file.read_text());assert row['guide']==g['guide'] and row['contract_sha256']==contract and row['sample_sha256']==sample_sha and sha(null)==row['null_SHA256']
            else:
                assert i>=106 and not null.exists()
                start=time.monotonic();seed=int.from_bytes(hashlib.sha256((str(cfg['permutation_seed'])+'|'+g['guide']).encode()).digest()[:8],'little')
                stat,values=test_pair(pc[g['rows']],pc[ref],seed,cfg['permutations'],official)
                tmp=null.with_name(null.name+'.receiving')
                with tmp.open('wb') as f:np.savez_compressed(f,null_energy_unbiased=values)
                tmp.replace(null)
                row=dict(guide=g['guide'],role=g['role'],target_gene_name=str(catalogue.loc[g['guide'],'target_gene_name']),target_gene_id=str(catalogue.loc[g['guide'],'target_gene_id']),
                    **stat,seconds=time.monotonic()-start,contract_sha256=contract,sample_sha256=sample_sha,query_barcodes=bc[g['rows']].tolist(),reference_barcodes=bc[ref].tolist(),
                    null_SHA256=sha(null),permutation_seed=seed,reference_seed=cfg['reference_seed'],median_RNA_UMI=float(np.median(totals[g['rows']])),
                    diagnostic_extension_commit=args.extension_commit,screening_reliable=False)
                atomic(file,row)
            results.append(row);progress('descriptive_guide_tests',completed=len(results),total=785,last_guide=g['guide'])
        pvalues=np.asarray([r['pvalue'] for r in results]);reject,adjusted,_,_=multipletests(pvalues,alpha=.05,method='fdr_bh')
        small=[]
        for row,padj,flag in zip(results,adjusted,reject):
            item={k:v for k,v in row.items() if k not in ['query_barcodes','reference_barcodes']};item.update(pvalue_adj=float(padj),significant_adj=bool(flag),screening_reliable=False);small.append(item)
        frame=pd.DataFrame(small);frame.to_csv(root/'tested_guides_diagnostic.csv',index=False)
        all_table=catalogue.reset_index().merge(frame,on='guide',how='left',suffixes=('_design',''));all_table['n_QC_records']=all_table.guide.map(counts).fillna(0).astype(int)
        all_table['assessment_status']=np.where(all_table.n_QC_records.eq(0),'no_QC_records',np.where(all_table.n_QC_records.lt(10),'below10_not_tested','descriptive_test_NTC_review_required'))
        all_table.to_csv(root/'all_design_guides_diagnostic.csv',index=False)
        summary=dict(task_id=task['task_id'],primary_statistic_commit=cfg['commit'],extension_commit=args.extension_commit,contract_sha256=contract,records=71664,tested_guides=785,all_design_guides=1196,
            calibration=diagnostic,reference_sensitivity=json.loads((root/'reference_sensitivity/summary.json').read_text()),screening_reliable=False,reliable_candidate_list_produced=False,
            targeting_nominal_significant=int((frame.role.isin(['targeting','positive_control'])&frame.significant_adj).sum()),NTC_nominal_significant=int((frame.role.eq('non_targeting_control')&frame.significant_adj).sum()),
            alpha=.05,permutations=cfg['permutations'],multiple_testing='BH across785 groups',estimator='equal-N unbiased energy U, same ranks as released pertpy V statistic',
            inference='descriptive RNA-state association in company merged records; NTC background differences preclude attribution to intended targets',functional_KO_efficiency_assessed=False)
        atomic(root/'diagnostic_summary.json',summary)
        atomic(root/'diagnostic_result_manifest.json',{n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in ['tested_guides_diagnostic.csv','all_design_guides_diagnostic.csv','diagnostic_summary.json','NTC_diagnostic.json']})
        progress('descriptive_tests_complete_pending_interpretation',completed=785,total=785)
    except BaseException as e:progress('failed_pending_recovery',error=repr(e));raise
