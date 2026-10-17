"""Correct undefined control-z testing domain without repeating unaffected work."""
from pathlib import Path
from datetime import datetime
import argparse,json,os,subprocess,time,shutil
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from analyze_expression_structure import sha256,atomic_json
from run_chip01_literature import verify_stage
from replogle2022_statistics import reference_de_domain


def verified_copy(source,target):
    if target.exists():raise FileExistsError(target)
    before=sha256(source);shutil.copyfile(source,target)
    if sha256(source)!=before or sha256(target)!=before:raise ValueError('Source/destination copy hash mismatch')
    return dict(bytes=target.stat().st_size,sha256=before)


def run(args):
    if args.output.exists():raise FileExistsError('New corrected output required')
    parent=verify_stage(args.source_effects)
    identity=json.loads((args.source_effects/'summary.json').read_text())
    genes=pd.read_csv(args.source_effects/'gene_statistics.csv.gz',index_col=0)
    domain=genes.relative_z_defined.to_numpy(dtype=bool)
    if not domain.any():raise ValueError('Empty control-relative z domain')
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**extra):atomic_json(args.output/'status.json',dict(stage=stage,pid=os.getpid(),time=datetime.now().astimezone().isoformat(),code_commit=commit,**extra))
    try:
        copies={}
        for name in ['cell_policy.csv.gz','energy_representation.npz','energy_results.csv','energy_progress.csv','direct_target_masking.csv']:
            copies[name]=verified_copy(args.source_effects/name,args.output/name)
        parameters=json.loads((args.source_effects/'frozen_parameters.json').read_text())
        parameters['de_gene_domain']='finite positive core-control SD; undefined relative-z genes stored but not tested'
        parameters['de_bh_family']='all and only defined control-z RNA genes per perturbation'
        atomic_json(args.output/'frozen_parameters.json',parameters)
        summaries=pd.read_csv(args.source_effects/'de_summary.csv');tops=[]
        for i,row in summaries.iterrows():
            table=pd.read_csv(args.source_effects/row.file)
            if not np.array_equal(table.gene_id.astype(str).to_numpy(),genes.index.astype(str).to_numpy()):raise ValueError('Gene identity/order mismatch')
            if not np.array_equal(table.response_z_defined.to_numpy(dtype=bool),domain):raise ValueError('Domain identity mismatch')
            table['de_gene_testable']=domain;table['ad_tested']=domain
            table['ad_pvalue'],table['ad_bh_q']=reference_de_domain(table.ad_pvalue.to_numpy(),domain)
            table.loc[~domain,'ad_statistic']=np.nan;table['ad_test_informative']=domain
            mw_tested=bool(row.mw_tested);table['mw_tested']=domain&mw_tested
            if mw_tested:table['mw_pvalue'],table['mw_bh_q']=reference_de_domain(table.mw_pvalue.to_numpy(),domain)
            else:table['mw_pvalue']=np.nan;table['mw_bh_q']=np.nan
            table.to_csv(args.output/row.file,index=False)
            summaries.loc[i,'tested_genes']=int(domain.sum())
            summaries.loc[i,'ad_de_genes']=int((table.ad_bh_q<.05).sum())
            summaries.loc[i,'mw_de_genes']=int((table.mw_bh_q<.05).sum()) if mw_tested else np.nan
            hits=table[table.ad_bh_q<.05].nlargest(50,'ad_statistic').copy()
            if len(hits):hits.insert(0,'target',row.target);tops.append(hits)
            if (i+1)%8==0 or i+1==len(summaries):state('correct_defined_z_de_domain',completed_targets=i+1,total_targets=len(summaries),tested_genes=int(domain.sum()),stored_genes=len(genes))
        summaries.to_csv(args.output/'de_summary.csv',index=False)
        if tops:pd.concat(tops,ignore_index=True).to_csv(args.output/'top_ad_genes.csv.gz',index=False)
        else:pd.DataFrame(columns=['target','gene_id','gene_name','ad_pvalue','ad_bh_q']).to_csv(args.output/'top_ad_genes.csv.gz',index=False)
        genes['de_gene_testable']=domain;genes.to_csv(args.output/'gene_statistics.csv.gz')
        with np.load(args.source_effects/'response_profiles.npz') as original:
            profiles={name:original[name] for name in original.files}
        if not np.array_equal(profiles['z_defined'],domain):raise ValueError('Profile domain identity mismatch')
        profiles['ad_statistics'][:,~domain]=-np.inf
        np.savez_compressed(args.output/'response_profiles.npz',**profiles)
        atomic_json(args.output/'de_domain_correction.json',dict(source_effects=str(args.source_effects),source_effects_commit=parent['code_commit'],
            source_effects_manifest_sha256=sha256(args.source_effects/'output_manifest.json'),corrected_commit=commit,
            tested_genes=int(domain.sum()),stored_genes=len(genes),excluded_undefined_control_z=int((~domain).sum()),
            energy_recomputed=False,core_controls_recomputed=False,verified_unchanged_copies=copies,
            reason='Control-relative z requires positive finite control SD; rank equivalence does not apply to zero SD',
            asymptotic_pvalue_null_calibration='Not established for extremely rare/tied expression; reference-calculation equivalence is not empirical false-positive validation',source_deleted=False))
        identity['tested_rna_features']=int(domain.sum());identity['de_gene_domain']='defined_core_control_relative_z'
        identity['source_effects_commit']=parent['code_commit'];identity['ad_targets_with_de']=int((summaries.ad_de_genes>0).sum())
        atomic_json(args.output/'summary.json',identity)
        atomic_json(args.output/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in args.output.iterdir() if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json']})
        state('complete',elapsed_seconds=time.monotonic()-start,source_modified=False,stage_completed='effects',source_effects_commit=parent['code_commit'])
    except Exception as exc:state('failed',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-effects',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    with threadpool_limits(2):run(args)
