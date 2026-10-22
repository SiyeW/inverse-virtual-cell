"""Staged CPU literature-aligned redo, immutable source and distinct outputs."""
from pathlib import Path
from datetime import datetime
import argparse,csv,hashlib,json,os,subprocess,time,importlib.metadata as md
import numpy as np
import pandas as pd
import anndata as ad
import h5py,psutil
from scipy import sparse
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
from analyze_expression_structure import blocks,sha256,atomic_json
from perturbation_statistics import bh_adjust,differential_expression_block
from replogle2022_statistics import (ADCalibration,ad_statistics,energy_permutation_full_control,reference_de_domain,
                                    strong_perturbation_mask,response_feature_union,paper_hdbscan)


def read_input(source,expected):
    if sha256(source)!=expected:raise ValueError('Source SHA256 mismatch')
    a=ad.read_h5ad(source,backed='r');obs=a.obs.copy();var=a.var.copy();shape=a.shape;a.file.close()
    selected=obs.qc_pass_stage1.astype(bool)&obs.guide_status.astype(str).eq('single_guide')
    if not np.array_equal(selected,obs.main_analysis_stage1.astype(bool)):raise ValueError('Existing candidate policy mismatch')
    rows=np.flatnonzero(selected.to_numpy())
    with h5py.File(source,'r') as f:
        if f['X'].attrs['encoding-type']!='csr_matrix':raise ValueError('Sparse CSR source required')
        counts=sparse.vstack([x for _,x in blocks(f,rows,len(var))],format='csr')
        nnz=len(f['X/data'])
    meta=obs.iloc[rows].copy();meta['source_row']=rows
    if not np.array_equal(np.asarray(counts.sum(1)).ravel(),meta.total_counts.to_numpy()):raise ValueError('Input count invariants differ')
    if (counts.data<0).any() or np.any(counts.data!=np.floor(counts.data)):raise ValueError('Expected original integer counts')
    return counts,meta,var,dict(source_shape=list(shape),source_nnz=nnz,selected_cells=len(rows),sha256=expected,source_modified=False)


def validate_calibration(calibration,validation):
    table=pd.read_csv(validation)
    observed=calibration.pvalues(table.statistic)
    error=np.max(np.abs(observed-table.pvalue.to_numpy())/np.maximum(table.pvalue.to_numpy(),1e-300))
    if error>.001:raise ValueError(f'R AD off-grid relative error {error} exceeds 0.1% gate')
    return dict(reference_points=len(table),maximum_relative_error=float(error),calibration_validated=True)


def validate_core_policy_identity(policy,original):
    """Compare identities, not pandas string-storage implementation details."""
    if policy.index.has_duplicates or original.index.has_duplicates or policy.index.isna().any() or original.index.isna().any():
        raise ValueError('Nonunique/missing source barcode')
    if not np.array_equal(policy.index.astype(str).to_numpy(),original.index.astype(str).to_numpy()):
        raise ValueError('Core/source barcode values or order mismatch')
    for column in ['source_row','total_counts','is_control']:
        if not np.array_equal(policy[column].to_numpy(),original[column].to_numpy()):
            raise ValueError('Core/source metadata mismatch: '+column)
    for column in ['guide_id','perturbed_gene']:
        if not np.array_equal(policy[column].astype(str).to_numpy(),original[column].astype(str).to_numpy()):
            raise ValueError('Core/source assignment mismatch: '+column)


def gene_symbols(var):
    """Categorical gene annotations must leave category storage before fillna."""
    return var.gene_name.astype('string').fillna('').to_numpy(dtype=str)


def verify_stage(path):
    status=json.loads((path/'status.json').read_text())
    if status['stage']!='complete':raise ValueError('Dependent stage is not complete')
    manifest=json.loads((path/'output_manifest.json').read_text())
    for rel,e in manifest.items():
        p=path/rel
        if p.stat().st_size!=e['bytes'] or sha256(p)!=e['sha256']:raise ValueError('Dependent artifact differs: '+rel)
    return status


def core_controls(args,config,state):
    state('verify_calibration')
    calibration=ADCalibration(args.calibration)
    validation=validate_calibration(calibration,args.calibration_validation)
    atomic_json(args.output/'calibration_validation.json',validation)
    state('verify_and_read_sparse_source')
    counts,obs,var,identity=read_input(args.source,args.input_sha256)
    control=obs.is_control.to_numpy(dtype=bool)
    guides=obs.guide_id.astype(str).to_numpy()
    ctrl=counts[control]
    means=np.asarray(ctrl.mean(0)).ravel()
    features=np.flatnonzero(means>config['core_control_gene_mean_raw_umi_gt'])
    if not len(features):raise ValueError('No genes satisfy paper core-control expression threshold')
    sizes=pd.Series(guides[control]).value_counts()
    candidates=sorted(sizes[sizes>=config['core_control_min_cells_for_pairwise_ad']].index)
    if len(candidates)<2:raise ValueError('Too few mathematically testable core-control candidates')
    sorted_groups={g:np.sort(counts[control&(guides==g)][:,features].toarray().astype(float),axis=0) for g in candidates}
    total=len(candidates)*(len(candidates)-1)//2;completed=0;started=time.monotonic()
    sums={g:0 for g in candidates};pair_counts={g:0 for g in candidates}
    with (args.output/'control_pairwise_ad.csv').open('w',encoding='utf8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['guide_a','guide_b','cells_a','cells_b','tested_genes','ad_bh_de_genes']);writer.writeheader()
        for i,left in enumerate(candidates):
            for right in candidates[i+1:]:
                statistics=ad_statistics(sorted_groups[left],sorted_groups[right],already_sorted=True)
                p=calibration.pvalues(statistics);q=bh_adjust(p);de=int((q<config['gene_bh_q']).sum())
                writer.writerow(dict(guide_a=left,guide_b=right,cells_a=len(sorted_groups[left]),cells_b=len(sorted_groups[right]),tested_genes=len(features),ad_bh_de_genes=de))
                sums[left]+=de;sums[right]+=de;pair_counts[left]+=1;pair_counts[right]+=1;completed+=1
                if completed%16==0 or completed==total:
                    stream.flush()
                    state('core_control_pairwise_ad',completed_pairs=completed,total_pairs=total,
                          tested_genes=len(features),elapsed_seconds=time.monotonic()-started,
                          seconds_per_pair=(time.monotonic()-started)/completed)
    records=[]
    for guide,size in sizes.items():
        testable=guide in pair_counts
        avg=sums[guide]/pair_counts[guide] if testable else None
        keep=testable and avg<config['core_control_mean_ad_de_genes_lt']
        records.append(dict(guide_id=guide,cells=int(size),pairwise_ad_testable=testable,mean_ad_de_genes=avg,core_control=keep,
                            reason='below_paper_core_threshold' if keep else ('above_paper_core_threshold' if testable else 'insufficient_cells_for_pairwise_AD')))
    summary=pd.DataFrame(records).sort_values('guide_id')
    summary.to_csv(args.output/'core_control_guides.csv',index=False)
    core_guides=set(summary.loc[summary.core_control,'guide_id'])
    obs['core_control']=control&np.isin(guides,list(core_guides))
    if obs.core_control.sum()<25:raise ValueError('Fewer than 25 core control cells; inspect results, do not change threshold automatically')
    obs['literature_analysis_cell']=~control|obs.core_control.to_numpy()
    obs.to_csv(args.output/'cell_policy.csv.gz')
    var.iloc[features].to_csv(args.output/'core_control_tested_genes.csv',index=False)
    identity.update(core_control_guides=len(core_guides),original_control_guides=len(sizes),core_control_cells=int(obs.core_control.sum()),
                    candidate_guides=len(candidates),testable_pairs=total,core_threshold=config['core_control_mean_ad_de_genes_lt'],
                    normalization_target_umi=float(np.median(obs.loc[obs.core_control,'total_counts'])),
                    calibration_sha256=sha256(args.calibration),calibration_validation=validation)
    atomic_json(args.output/'core_summary.json',identity)


def prepare_and_effects(args,config,state):
    core_state=verify_stage(args.core)
    calibration=ADCalibration(args.calibration)
    validation=validate_calibration(calibration,args.calibration_validation)
    state('verify_and_read_sparse_source')
    counts,original,var,identity=read_input(args.source,args.input_sha256)
    policy=pd.read_csv(args.core/'cell_policy.csv.gz',index_col=0)
    validate_core_policy_identity(policy,original)
    keep=policy.literature_analysis_cell.to_numpy(dtype=bool)
    obs=policy.loc[keep].copy();counts=counts[keep]
    controls=obs.core_control.to_numpy(dtype=bool);labels=obs.perturbed_gene.astype(str).to_numpy()
    if args.smoke:
        sizes=pd.Series(labels[~controls]).value_counts()
        picked=sorted(sizes[(sizes>=25)&(sizes<=100)].index)[:3]
        if len(picked)!=3:raise ValueError('Insufficient three bounded real-data smoke targets')
        select=controls|np.isin(labels,picked)
        obs=obs.iloc[np.flatnonzero(select)].copy();counts=counts[select];controls=controls[select];labels=labels[select]
        feature_mean=np.asarray(counts.mean(0)).ravel()
        subset=np.argsort(feature_mean,kind='stable')[-64:]
        counts=counts[:,subset];var=var.iloc[subset].copy()
    if not np.array_equal(controls,obs.is_control.to_numpy(dtype=bool)):raise ValueError('Non-core NTC leaked into new primary population')
    reference_total=float(np.median(obs.loc[controls,'total_counts']))
    totals=obs.total_counts.to_numpy(dtype=float)
    norm=(sparse.diags(reference_total/totals)@counts.astype(float)).tocsr()
    raw_mean=np.asarray(counts.mean(0)).ravel()
    raw_variance=np.maximum(0,np.asarray(counts.astype(float).power(2).mean(0)).ravel()-raw_mean**2)
    ctrl=norm[controls].tocsc();n_control=int(controls.sum())
    center=np.asarray(ctrl.mean(0)).ravel()
    scale=np.sqrt(np.maximum(0,(np.asarray(ctrl.power(2).mean(0)).ravel()-center**2)*n_control/(n_control-1)))
    defined=np.isfinite(scale)&(scale>1e-12)
    targets=sorted(set(labels[~controls]));indices=[np.flatnonzero((labels==t)&~controls) for t in targets]
    cells=np.array([len(i) for i in indices]);genes=len(var)
    ad_scores=np.empty((len(targets),genes));ad_p=np.empty_like(ad_scores)
    mw_p=np.ones_like(ad_scores);means=np.empty_like(ad_scores);fractions=np.empty_like(ad_scores)
    response=np.zeros_like(ad_scores);rank_effect=np.zeros_like(ad_scores)
    normalized=norm.tocsc();started=time.monotonic();chunks=(genes+127)//128
    for block,start in enumerate(range(0,genes,128)):
        end=min(start+128,genes)
        dense=normalized[:,start:end].toarray()
        control=dense[controls];sorted_control=np.sort(control,axis=0)
        for i,picked in enumerate(indices):
            x=dense[picked];means[i,start:end]=x.mean(0);fractions[i,start:end]=(x>0).mean(0)
            statistic=ad_statistics(np.sort(x,axis=0),sorted_control,already_sorted=True)
            ad_scores[i,start:end]=statistic;ad_p[i,start:end]=calibration.pvalues(statistic)
            if len(x)>=config['de_mw_min_cells']:
                mw=differential_expression_block(x,control)
                mw_p[i,start:end]=mw['pvalue'];rank_effect[i,start:end]=mw['rank_biserial']
        state('literature_gene_statistics',completed_gene_blocks=block+1,total_gene_blocks=chunks,
              completed_genes=end,total_genes=genes,targets=len(targets),control_cells=n_control,
              elapsed_seconds=time.monotonic()-started,seconds_per_gene_block=(time.monotonic()-started)/(block+1))
    response[:,defined]=(means[:,defined]-center[defined])/scale[defined]
    names=gene_symbols(var)
    direct_mask=[]
    for i,target in enumerate(targets):
        matched=np.flatnonzero(names==target)
        masked=len(matched)==1 and defined[matched[0]]
        if masked:response[i,matched[0]]=0
        direct_mask.append(dict(target=target,exact_target_rna_feature_masked=bool(masked),ambiguous_matches=len(matched)))
    pd.DataFrame(direct_mask).to_csv(args.output/'direct_target_masking.csv',index=False)
    summaries=[];top=[];summary_started=time.monotonic()
    control_fraction=np.asarray((ctrl>0).mean(0)).ravel()
    for i,target in enumerate(targets):
        ad_p[i],aq=reference_de_domain(ad_p[i],defined)
        if cells[i]>=config['de_mw_min_cells']:mw_p[i],mq=reference_de_domain(mw_p[i],defined)
        else:mw_p[i]=np.nan;mq=np.full(genes,np.nan)
        ad_scores[i,~defined]=-np.inf;defined_ad=np.isfinite(ad_scores[i])
        table=pd.DataFrame(dict(gene_id=var.index.astype(str),gene_name=names,
            ad_statistic=np.where(defined_ad,ad_scores[i],np.nan),ad_test_informative=defined_ad,
            de_gene_testable=defined,ad_tested=defined,
            ad_pvalue=ad_p[i],ad_bh_q=aq,mw_tested=defined&(cells[i]>=config['de_mw_min_cells']),
            mw_pvalue=mw_p[i],mw_bh_q=mq,rank_biserial=rank_effect[i],
            mean_normalized=means[i],control_mean_normalized=center,
            log2_fold_change=np.log2((means[i]+1)/(center+1)),fraction_detected=fractions[i],
            control_fraction_detected=control_fraction,mean_response_z=response[i],response_z_defined=defined))
        file=f'target_{i:03d}_de.csv.gz';table.to_csv(args.output/file,index=False)
        summaries.append(dict(target=target,file=file,cells=int(cells[i]),control_cells=n_control,tested_genes=int(defined.sum()),
                              ad_de_genes=int((aq<config['gene_bh_q']).sum()),mw_tested=bool(cells[i]>=config['de_mw_min_cells']),
                              mw_de_genes=int((mq<config['gene_bh_q']).sum()) if cells[i]>=config['de_mw_min_cells'] else None,
                              low_coverage_for_clustering=bool(cells[i]<25)))
        hits=table[table.ad_bh_q<config['gene_bh_q']].sort_values('ad_statistic',ascending=False).head(50).copy()
        if len(hits):hits.insert(0,'target',target);top.append(hits)
    pd.DataFrame(summaries).to_csv(args.output/'de_summary.csv',index=False)
    if top:pd.concat(top,ignore_index=True).to_csv(args.output/'top_ad_genes.csv.gz',index=False)
    else:pd.DataFrame(columns=['target','gene_id','gene_name','ad_pvalue','ad_bh_q']).to_csv(args.output/'top_ad_genes.csv.gz',index=False)
    gene_statistics=var.copy();gene_statistics['raw_mean']=raw_mean;gene_statistics['raw_variance']=raw_variance
    gene_statistics['core_control_normalized_mean']=center;gene_statistics['core_control_normalized_sd']=scale
    gene_statistics['relative_z_defined']=defined
    gene_statistics['de_gene_testable']=defined
    gene_statistics.to_csv(args.output/'gene_statistics.csv.gz')
    np.savez_compressed(args.output/'response_profiles.npz',targets=np.array(targets),values=response,
                        ad_statistics=ad_scores,raw_mean=raw_mean,raw_variance=raw_variance,z_defined=defined)
    obs.to_csv(args.output/'cell_policy.csv.gz')
    state('literature_energy_representation')
    features=np.flatnonzero((raw_mean>config['energy_gene_mean_raw_umi_gt'])&defined)
    if min(len(features),len(obs))<=config['energy_pcs']:raise ValueError('Insufficient features for reference energy PCA')
    if len(obs)*len(features)*8>4*1024**3:raise MemoryError('Selected energy feature matrix exceeds explicit 4GiB bound')
    z=norm[:,features].toarray();z=(z-center[features])/scale[features]
    np.minimum(z,config['energy_upper_z_clip'],out=z)
    pca=PCA(n_components=config['energy_pcs'],svd_solver='randomized',random_state=config['seed'])
    scores=pca.fit_transform(z);del z
    reference=np.flatnonzero(controls)
    if len(reference)>config['energy_max_controls']:
        reference=np.sort(np.random.default_rng(config['seed']).choice(reference,config['energy_max_controls'],replace=False))
    np.savez_compressed(args.output/'energy_representation.npz',coordinates=scores,feature_indices=features,
                        source_rows=obs.source_row.to_numpy(),reference_rows=reference,explained_variance_ratio=pca.explained_variance_ratio_)
    energies=[];energy_started=time.monotonic()
    for i,target in enumerate(targets):
        effect,p=energy_permutation_full_control(scores[indices[i]],scores[reference],
                    permutations=config['energy_label_permutations'],seed=config['seed']+i)
        energies.append(dict(target=target,cells=int(cells[i]),control_cells=len(reference),energy=effect,permutation_p=p,
                             permutations=config['energy_label_permutations'],low_coverage=bool(cells[i]<25)))
        pd.DataFrame(energies).to_csv(args.output/'energy_progress.csv',index=False)
        state('literature_energy_permutations',completed_targets=i+1,total_targets=len(targets),target=target,
              control_cells=len(reference),elapsed_seconds=time.monotonic()-energy_started,
              seconds_per_target=(time.monotonic()-energy_started)/(i+1))
    energy=pd.DataFrame(energies);energy['bh_q']=bh_adjust(energy.permutation_p.to_numpy())
    energy.to_csv(args.output/'energy_results.csv',index=False)
    identity.update(analysis_cells=len(obs),core_control_cells=n_control,targets=len(targets),rna_features=genes,
                    tested_rna_features=int(defined.sum()),de_gene_domain='defined_core_control_relative_z',
                    energy_features=len(features),energy_control_cells=len(reference),
                    ad_targets_with_de=int(sum(r['ad_de_genes']>0 for r in summaries)),
                    energy_raw_p_lt_005=int((energy.permutation_p<.05).sum()),energy_bh_q_lt_005=int((energy.bh_q<.05).sum()),
                    reference_core_commit=core_state['code_commit'],calibration_validation=validation,
                    scientific_result=not args.smoke,engineering_smoke=bool(args.smoke))
    atomic_json(args.output/'summary.json',identity)


def run(args):
    if args.output.exists():raise FileExistsError('Distinct run output required; inspect before recovery')
    if psutil.virtual_memory().available<12*1024**3:raise MemoryError('Execution host needs 12GiB free; do not run full counts locally')
    config=json.loads(args.config.read_text())
    if args.smoke:
        if args.stage!='effects':raise ValueError('Smoke supported only for effects')
        config['engineering_smoke']=True;config['energy_label_permutations']=99
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**kw):atomic_json(args.output/'status.json',dict(stage=stage,pid=os.getpid(),time=datetime.now().astimezone().isoformat(),code_commit=commit,**kw))
    try:
        atomic_json(args.output/'frozen_parameters.json',config)
        with threadpool_limits(config['threads']):
            if args.stage=='core':core_controls(args,config,state)
            elif args.stage=='effects':prepare_and_effects(args,config,state)
            else:raise ValueError('Unsupported stage')
        atomic_json(args.output/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in args.output.iterdir() if p.is_file() and p.name not in ['status.json','owner.lock','output_manifest.json']})
        state('complete',elapsed_seconds=time.monotonic()-start,stage_completed=args.stage,
              versions={k:md.version(k) for k in ['numpy','scipy','numba','anndata','scikit-learn']},source_modified=False)
    except Exception as error:state('failed',error=repr(error));raise
    finally:lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['core','effects'],required=True)
    for name in ['source','output','config','calibration','calibration-validation']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--input-sha256',required=True)
    parser.add_argument('--core',type=Path)
    parser.add_argument('--smoke',action='store_true')
    run(parser.parse_args())
