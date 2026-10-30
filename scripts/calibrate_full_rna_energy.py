"""Calibrate complete-population RNA energy tests on real NTC guide groups.

No targeting guide is selected here. Independent label permutations are
compared with a prespecified RNA-depth-blocked variant; the latter asks about
differences beyond coarse measured-depth composition, not absolute RNA amount.
This is an adaptation of energy permutation testing, not the SCEPTRE pipeline.
"""
import argparse,atexit,hashlib,json,os,socket,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from rna_energy_permutation import subset_statistics,permutation_samples,depth_strata,permutation_pvalue

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def atomic(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.receiving')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2))
    for attempt in range(5):
        try:tmp.replace(path);return
        except PermissionError:
            if attempt==4:raise
            time.sleep(.15*(attempt+1))

def run(cfg):
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    op=Path(cfg['native']);out=Path(cfg['output']);out.mkdir(exist_ok=True)
    lock=op/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    st=Path('/proc/self/stat').read_text();ticks=int(st[st.rfind(')')+2:].split()[19])
    os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,commit=cfg['commit'])).encode());os.close(fd)
    atexit.register(lambda:lock.unlink(missing_ok=True))
    started=time.monotonic()
    def progress(stage,**kw):atomic(op/'status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),elapsed=time.monotonic()-started,**kw))
    try:
        progress('verify_small_accepted_inputs')
        for p,h in cfg.get('source_sha256',{}).items():
            if sha(p)!=h:raise ValueError('Fixed scientific source differs: '+p)
        for p,h in cfg['input_sha256'].items():
            if sha(p)!=h:raise ValueError('Input differs: '+p)
        with np.load(cfg['pc_input'],allow_pickle=False) as z:
            pc=z['PCA'];labels=z['guide'];roles=z['role'];barcodes=z['cell_barcode'];totals=z['total_RNA'];rna_ids=z['RNA_ids']
        assert pc.shape==(71664,40) and len(set(barcodes))==71664 and len(np.unique(labels))==785
        assert hashlib.sha256(json.dumps(rna_ids.tolist()).encode()).hexdigest()==cfg['rna_inventory_sha256']
        nc=np.flatnonzero(roles=='non_targeting_control');assert len(nc)==19342
        x=np.asarray(pc[nc],dtype=np.float64);guide=labels[nc];depth=totals[nc]
        del pc
        old=pd.read_csv(cfg['accepted_effect_table']).set_index('guide',verify_integrity=True)
        progress('NTC_complete_pair_geometry',records=len(nc),RAM_distance_bytes=int(len(nc)**2*8))
        d=cdist(x,x);np.fill_diagonal(d,0.);row=d.sum(axis=1);total=float(row.sum())
        empty_cross=np.empty((0,len(nc)),dtype=float);empty_query=np.empty((0,0),dtype=float)
        blocks=depth_strata(depth,10);groups=sorted(np.unique(guide));assert len(groups)==106
        contract=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()
        atomic(out/'input_proof.json',dict(config_contract=contract,records=71664,NTC_records=19342,NTC_guides=106,
            full_query=True,full_NTC=True,real_reference_excludes_whole_query_guide=True,RNA_scope='frozen HVG2000/PCA40',
            test_statistic='released pertpy energy V; standard unbiased U is reported separately',
            conditional_variant='permutations within label-blind deciles of measured total RNA; no group downsampling',
            inference_scope='company_merged_records_diagnostic; not capture-resolved biological effects'))
        results=[];dest=out/'NTC_guides';dest.mkdir(exist_ok=True)
        for number,g in enumerate(groups):
            file=dest/(hashlib.sha256(str(g).encode()).hexdigest()[:20]+'.json')
            if file.exists():
                result=json.loads(file.read_text());assert result['guide']==str(g) and result['contract']==contract
            else:
                selected=np.flatnonzero(guide==g).astype('int32');n=len(selected);m=len(nc)-n
                obs=subset_statistics(selected[None,:],d,empty_cross,empty_query,row,total)[0]
                expected=float(old.loc[str(g),'energy_unbiased'])
                if not np.isclose(obs[1],expected,atol=1e-8,rtol=1e-7):raise ValueError('Accepted full-reference energy bridge failed: '+str(g))
                result=dict(guide=str(g),n_query=n,n_reference=m,energy_V=float(obs[0]),energy_U=float(obs[1]),accepted_energy_U=expected,contract=contract)
                for method,strata in [('independent',None),('RNA_depth_deciles',blocks)]:
                    seed=int.from_bytes(hashlib.sha256((str(cfg['seed'])+'|'+method+'|'+str(g)).encode()).digest()[:8],'little')
                    phase=time.monotonic();draws=permutation_samples(len(nc),selected,cfg['permutations'],seed,strata)
                    null=subset_statistics(draws,d,empty_cross,empty_query,row,total)
                    p,extreme=permutation_pvalue(obs[0],null[:,0])
                    result[method]=dict(pvalue=p,extreme=extreme,permutations=cfg['permutations'],seconds=time.monotonic()-phase,
                        null_U_mean=float(null[:,1].mean()),null_U_SD=float(null[:,1].std(ddof=1)))
                atomic(file,result)
            results.append(result)
            progress('NTC_full_reference_calibration',completed=len(results),total=106,last_guide=str(g),last_n=int(result['n_query']),
                     last_guide_seconds=sum(result[k]['seconds'] for k in ['independent','RNA_depth_deciles']))
        rows=[]
        for r in results:
            for method in ['independent','RNA_depth_deciles']:
                rows.append({k:r[k] for k in ['guide','n_query','n_reference','energy_V','energy_U','accepted_energy_U']}|dict(method=method,**r[method]))
        frame=pd.DataFrame(rows);frame.to_csv(out/'NTC_calibration.csv',index=False)
        summary={}
        for method,g in frame.groupby('method'):
            p=g.pvalue.to_numpy();flag=int((p<.05).sum())
            summary[method]=dict(tests=len(p),raw_P_lt05=flag,fraction=flag/len(p),p_quantiles=np.quantile(p,[0,.1,.25,.5,.75,.9,1]).tolist(),
                                 meets_prespecified_review_guard=bool(flag/len(p)<=cfg['NTC_review_guard_fraction']))
        atomic(out/'calibration_summary.json',dict(methods=summary,review_guard_fraction=cfg['NTC_review_guard_fraction'],
            guard_is_diagnostic_not_proof_of_biological_FDR=True,targeting_tests_run=False,candidate_list_produced=False,
            bridge_max_abs_error=float(np.max(np.abs(frame.energy_U-frame.accepted_energy_U))),config_contract=contract))
        progress('NTC_calibration_complete_pending_method_judgment',completed=106,total=106)
    except BaseException as e:
        progress('failed_pending_recovery',error=repr(e));raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    run(json.loads(args.config.read_text()))
