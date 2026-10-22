"""Control-aware RNA QC and ordered target cards from immutable sparse counts."""
from pathlib import Path
from datetime import datetime
import argparse,hashlib,importlib.metadata as md,json,os,subprocess,time,traceback
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
from analyze_expression_structure import sha256
from copy_verified_delivery import atomic_json
from run_chip01_literature import read_input,validate_core_policy_identity
from chip01_qc_supplement import library_normalize,distribution_summary,permutation_rank_diagnostic
from run_chip01_qc_supplement import negative_control_checks
from chip01_report_alignment import (TYPES,annotation_overlay,reference_masks,effect_profiles,
    expression_effect,leave_guide_out_reference,stable_rng,target_support,qc_flow)
from plot_chip01_report_alignment import (COLORS,NAMES,technical_figures,reference_figure,
    effect_background_figures,target_atlas)


def read_prefix_smoke(source,all_obs,var,prefix_rows):
    """A bounded contiguous CSR prefix; formal run still hashes/reads the full source."""
    stop=min(len(all_obs),int(prefix_rows))
    with h5py.File(source,'r') as f:
        x=f['X']
        if x.attrs['encoding-type']!='csr_matrix':raise ValueError('Bounded smoke requires sparse CSR source')
        pointer=x['indptr'][:stop+1];end=int(pointer[-1])
        matrix=sparse.csr_matrix((x['data'][:end],x['indices'][:end],pointer),shape=(stop,len(var)))
    if not np.array_equal(np.asarray(matrix.sum(1)).ravel(),all_obs.total_counts.to_numpy()[:stop]):
        raise ValueError('Smoke source count totals differ')
    keep=all_obs.main_analysis_stage1.astype(bool).to_numpy()[:stop]
    obs=all_obs.iloc[np.flatnonzero(keep)].copy();obs['source_row']=np.flatnonzero(keep)
    return matrix[keep],obs,dict(source_shape=[len(all_obs),len(var)],source_prefix_rows=stop,
        source_prefix_nnz=end,selected_cells=len(obs),source_full_hash_recomputed=False,source_modified=False)


def role_summary(overlay,records,output):
    rows=[]
    for kind in TYPES:
        part=overlay.loc[overlay.source_control_type.eq(kind)]
        values=part.single_guide_after_stage1_qc.to_numpy()
        rows.append(dict(source_control_type=kind,library_guides=len(part),cells=int(values.sum()),
                         eligible_guides=int(records.source_control_type.eq(kind).sum()),
                         **{k:v for k,v in distribution_summary(values,[1,10,20]).items() if k!='n'}))
    pd.DataFrame(rows).to_csv(output/'tables/control_type_summary.csv',index=False)
    return rows


def example_means(norm,obs,overlay,reference,matched_cells,seed):
    guides=obs.guide_id.astype(str).to_numpy();examples=[];selections=[]
    for kind in TYPES[:2]:
        eligible=overlay.loc[overlay.source_control_type.eq(kind)].sort_values(['single_guide_after_stage1_qc','guide'],ascending=[False,True])
        chosen=None
        for row in eligible.itertuples():
            ix=np.flatnonzero(guides==row.guide)
            if len(ix)>=matched_cells:chosen=(row,ix);break
        if chosen is None:continue
        row,ix=chosen;center,ref=leave_guide_out_reference(norm,guides,reference,row.guide)
        picked=np.sort(stable_rng(seed,'example:'+row.guide).choice(ix,matched_cells,replace=False))
        examples.append(dict(label=f'{NAMES[kind]}: highest-coverage guide; n={matched_cells}',mean=np.asarray(norm[picked].mean(0)).ravel(),reference_mean=center,color=COLORS[kind]))
        selections.append(dict(group=row.guide,role=kind,cells=matched_cells,source_rows=obs.iloc[picked].source_row.astype(int).tolist(),reference_cells=len(ref),reference_overlap_cells=0))
    center=np.asarray(norm[reference].mean(0)).ravel()
    for target,part in overlay.loc[overlay.source_control_type.eq('positive_control')].groupby('target_gene_name',sort=True):
        ix=np.flatnonzero(np.isin(guides,part.guide))
        if len(ix)<matched_cells:continue
        picked=np.sort(stable_rng(seed,'example:'+target).choice(ix,matched_cells,replace=False))
        if np.intersect1d(picked,np.flatnonzero(reference)).size:raise ValueError('Example/reference overlap')
        examples.append(dict(label=f'Annotated positive: {target}; n={matched_cells}; expected RNA phenotype unverified',mean=np.asarray(norm[picked].mean(0)).ravel(),reference_mean=center,color=COLORS['positive_control']))
        selections.append(dict(group=str(target),role='positive_control',cells=matched_cells,source_rows=obs.iloc[picked].source_row.astype(int).tolist(),reference_cells=int(reference.sum()),reference_overlap_cells=0))
    return examples,selections


def spot_pure_null(counts,obs,genes,config,output):
    """One bounded rare-tie diagnostic, not a replacement genome-wide test."""
    from scipy.stats import mannwhitneyu
    selection=json.loads((output/'null_checks/null_selection.json').read_text())
    chosen=selection['n25_r01'];lookup={str(x):i for i,x in enumerate(obs.index)}
    treatment=np.array([lookup[x] for x in chosen['treatment_barcodes']]);control=np.array([lookup[x] for x in chosen['control_barcodes']])
    features=np.flatnonzero(genes.relative_z_defined.to_numpy(dtype=bool))
    nt=counts[treatment][:,features];nc=counts[control][:,features]
    candidates=np.flatnonzero((np.asarray((nt>0).sum(0)).ravel()==1)&(np.asarray((nc>0).sum(0)).ravel()==0))
    if not len(candidates):return dict(assessed=False,reason='No one-versus-zero diagnostic pair in the fixed n25 comparison')
    j=sorted(candidates,key=lambda k:str(genes.index[features[k]]))[0];feature=features[j]
    x=library_normalize(counts[treatment,feature],obs.iloc[treatment].total_counts,config['normalization_target_umi']).toarray().ravel()
    y=library_normalize(counts[control,feature],obs.iloc[control].total_counts,config['normalization_target_umi']).toarray().ravel()
    p=np.load(output/'null_checks/tables/null_pvalues.npz')
    row=list(p['comparisons']).index('n25_r01')
    return dict(assessed=True,comparison='n25_r01',gene_id=str(genes.index[feature]),gene_name=str(genes.iloc[feature].gene_name),treatment_detected=1,control_detected=0,
                source_cells=len(treatment),control_cells=len(control),stored_ad_pvalue=float(p['ad_pvalues'][row,j]),stored_mw_pvalue=float(p['mw_pvalues'][row,j]),
                independent_scipy_mw_pvalue=float(mannwhitneyu(x,y,alternative='two-sided',method='asymptotic').pvalue),
                randomization=permutation_rank_diagnostic(x,y,19999,config['seed']),single_positive_allocation_probability=len(x)/(len(x)+len(y)),
                interpretation='Pure-NTC diagnostic only; no genome-wide p-value correction or biological hit claim')


def run(args):
    config=json.loads(args.config.read_text());out=Path(config['output']);runtime=Path(config['runtime']);source=Path(config['source'])
    if out.parent.resolve()!=Path(config['allowed_output_parent']).resolve() or out.exists() or runtime.exists():raise ValueError('Distinct approved output/runtime required; no automatic rerun')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip()
    if commit!=config['expected_commit']:raise ValueError('Execution commit differs from pushed fixed version')
    out.mkdir();runtime.mkdir(parents=True)
    for rel in ['tables','figures','figures/targets','null_checks/tables','null_checks/figures','null_runtime']:(out/rel).mkdir(exist_ok=True,parents=True)
    owner=dict(pid=os.getpid(),uid=os.getuid(),code_commit=commit,started_at=datetime.now().astimezone().isoformat())
    lock=out/'owner.lock';fd=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    with os.fdopen(fd,'w') as stream:json.dump(owner,stream)
    started=time.monotonic();before=source.stat()
    def state(stage,**more):
        value=dict(stage=stage,time=datetime.now().astimezone().isoformat(),elapsed_seconds=time.monotonic()-started,**owner,**more)
        atomic_json(runtime/'status.json',value);print(json.dumps(value),flush=True)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
        state('verify_small_inputs')
        inputs={}
        for key in ['coverage','core_policy','gene_statistics','annotation_audit','calibration']:
            p=Path(config[key]);inputs[key]=dict(path=str(p),bytes=p.stat().st_size,sha256=sha256(p))
            if inputs[key]['sha256']!=config['input_sha256'][key]:raise ValueError('Frozen input changed: '+key)
        atomic_json(out/'input_manifest.json',inputs)
        audit=json.loads(Path(config['annotation_audit']).read_text());coverage=pd.read_csv(config['coverage'])
        if 'guide_type' not in coverage:coverage['guide_type']=np.where(coverage.target_gene_name.eq('Non-Targeting'),'NTC','targeting')
        overlay=annotation_overlay(coverage,audit['annotations'])
        view=ad.read_h5ad(source,backed='r');all_obs=view.obs.copy();var=view.var.copy();view.file.close()
        if [len(all_obs),len(var)]!=config['expected_shape']:raise ValueError('Source dimensions differ')
        state('read_source_prefix' if args.smoke else 'verify_source_hash_and_read_sparse_counts')
        if args.smoke:counts,obs,identity=read_prefix_smoke(source,all_obs,var,config.get('smoke_prefix_rows',3000))
        else:counts,obs,var,identity=read_input(source,config['source_sha256'])
        policy=pd.read_csv(config['core_policy'],index_col=0).loc[obs.index]
        validate_core_policy_identity(policy,obs)
        old_core=policy.core_control.to_numpy(dtype=bool)
        cell_types,pure_ntc,reference=reference_masks(obs,overlay,old_core)
        if not args.smoke:
            if old_core.sum()!=config['expected_historical_core_cells'] or reference.sum()!=config['expected_pure_core_cells'] or pure_ntc.sum()!=config['expected_pure_ntc_cells']:
                raise ValueError('Frozen core/control populations differ')
        overlay['historical_core_guide']=overlay.guide.isin(obs.loc[old_core,'guide_id'].astype(str))
        overlay['pure_core_guide']=overlay.guide.isin(obs.loc[reference,'guide_id'].astype(str))
        overlay.to_csv(out/'tables/control_annotations.csv',index=False)
        genes=pd.read_csv(config['gene_statistics'],index_col=0)
        if not np.array_equal(genes.index.astype(str).to_numpy(),var.index.astype(str).to_numpy()):raise ValueError('RNA feature identity/order differs')
        features=np.flatnonzero((genes.raw_mean.to_numpy()>.5)&genes.relative_z_defined.to_numpy(dtype=bool))
        if len(features)!=config['expected_response_features']:raise ValueError('Frozen response feature family differs')
        if args.smoke:features=features[:64]
        names=genes.gene_name.astype('string').fillna('').to_numpy(dtype=str)
        state('reference_sensitivity')
        normalized=library_normalize(counts,obs.total_counts.to_numpy(),config['normalization_target_umi'])
        old_mean=np.asarray(normalized[old_core].mean(0)).ravel();pure_mean=np.asarray(normalized[reference].mean(0)).ravel()
        all_ntc_mean=np.asarray(normalized[pure_ntc].mean(0)).ravel()
        if not args.smoke and not np.allclose(old_mean,genes.core_control_normalized_mean.to_numpy(),rtol=1e-8,atol=1e-10):raise ValueError('Historical reference means cannot be independently reproduced')
        impact=pd.DataFrame(dict(gene_id=genes.index,gene_name=names,historical_mean=old_mean,pure_core_mean=pure_mean,all_pure_ntc_mean=all_ntc_mean,
            reference_log2_change=expression_effect(pure_mean,old_mean),all_ntc_vs_pure_core_log2_change=expression_effect(all_ntc_mean,pure_mean),
            frozen_response_feature=np.isin(np.arange(len(genes)),features)))
        impact.to_csv(out/'tables/reference_sensitivity.csv.gz',index=False)
        reference_figure(impact,int(old_core.sum()),int(reference.sum()),out)
        core=dict(historical_core_guides=int(overlay.historical_core_guide.sum()),historical_core_cells=int(old_core.sum()),
                  pure_core_guides=int(overlay.pure_core_guide.sum()),pure_core_cells=int(reference.sum()),
                  removed_intergenic_core_cells=int((old_core&~reference).sum()),all_pure_ntc_cells=int(pure_ntc.sum()),
                  median_abs_reference_log2_change=float(impact.reference_log2_change.abs().median()),
                  p99_abs_reference_log2_change=float(impact.reference_log2_change.abs().quantile(.99)),
                  max_abs_reference_log2_change=float(impact.reference_log2_change.abs().max()),normalization_target_umi=config['normalization_target_umi'],
                  all_ntc_vs_pure_core_median_abs_log2_change=float(impact.all_ntc_vs_pure_core_log2_change.abs().median()),
                  all_ntc_vs_pure_core_p99_abs_log2_change=float(impact.all_ntc_vs_pure_core_log2_change.abs().quantile(.99)),
                  descriptive_effect_reference='all source-annotated NTCs with whole-guide leave-out, not AD-selected core',
                  historical_outputs_modified=False)
        atomic_json(out/'reference_summary.json',core)
        selected=normalized[:,features].tocsr();del normalized
        minimum=2 if args.smoke else config['min_guide_cells'];matched=2 if args.smoke else config['matched_cells'];repeats=2 if args.smoke else config['matched_repeats']
        state('guide_effect_profiles',rna_features=len(features),matched_cells=matched,repeats=repeats)
        records,means,ref_means,effects,first,metrics,sampling=effect_profiles(selected,obs,overlay,pure_ntc,names[features],names,
             minimum,matched,repeats,config['seed'])
        records.to_csv(out/'tables/guide_effect_metadata.csv',index=False);metrics.to_csv(out/'tables/matched_effect_metrics.csv',index=False)
        atomic_json(out/'matched_cell_selections.json',sampling)
        np.savez_compressed(out/'tables/rna_effect_profiles.npz',means=means,reference_means=ref_means,log2fc=effects,matched_first=first,
            guides=records.guide.to_numpy(dtype=str),gene_ids=genes.index.to_numpy(dtype=str)[features],feature_indices=features,pseudocount=1.)
        roles=role_summary(overlay,records,out)
        state('technical_qc_and_background_figures')
        flow=qc_flow(all_obs);flow.to_csv(out/'tables/technical_cell_flow.csv',index=False);technical_figures(all_obs,overlay,flow,out)
        examples,selections=example_means(selected,obs,overlay,pure_ntc,matched,config['seed'])
        atomic_json(out/'example_cell_selections.json',selections)
        background=effect_background_figures(records,effects,first,metrics,examples,out)
        state('per_target_support_cards')
        support,pairs=target_support(selected,obs,overlay,records,effects,pure_ntc)
        support.to_csv(out/'tables/target_support.csv',index=False);pairs.to_csv(out/'tables/within_target_pairs.csv',index=False)
        cards=target_atlas(overlay,records,effects,support,pairs,out,max_cards=3 if args.smoke else None)
        state('pure_ntc_negative_control_checks')
        pure_obs=obs.copy();pure_obs['is_control']=pure_ntc
        null=negative_control_checks(counts,pure_obs,genes,reference,config,out/'null_checks',out/'null_runtime',state,args.smoke)
        spot=dict(assessed=False,reason='Engineering smoke') if args.smoke else spot_pure_null(counts,pure_obs,genes,config,out)
        atomic_json(out/'pure_ntc_spot_diagnostic.json',spot)
        after=source.stat()
        if (before.st_size,before.st_mtime_ns,before.st_dev,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_dev,after.st_ino):raise ValueError('Source identity changed during read-only analysis')
        for key,entry in inputs.items():
            if sha256(Path(entry['path']))!=entry['sha256']:raise ValueError('Selected input changed during analysis: '+key)
        summary=dict(run_id=config['run_id'],code_commit=commit,scientific_result=not args.smoke,engineering_smoke=bool(args.smoke),
            source_identity=identity,source_modified=False,core_reference=core,control_types=roles,response_features=len(features),
            eligible_guides=len(records),matched_cells=matched,matched_repeats=repeats,matched_effect_background=background.to_dict(orient='records'),
            targets=len(support),targets_assessable=int(support.assessment.eq('descriptive_assessable').sum()),
            targets_unassessable=int(support.assessment.ne('descriptive_assessable').sum()),target_cards=len(cards),
            negative_controls=null,spot_diagnostic=spot,cell_filtering_changed=False,pseudocount=1.,
            effect_reference='All source-annotated NTC cells, excluding intergenic controls; each NTC guide leaves all its own cells out of its reference; historical and pure core assessed separately for sensitivity',
            original_de_recalibrated=False,validated_functional_hits=False,new_clustering_run=False,dna_editing_efficiency='not_assessed',
            candidate_positive_control_expected_phenotype='not established for this experimental context',
            versions={k:md.version(k) for k in ['numpy','scipy','pandas','anndata','matplotlib']},settings=config,elapsed_seconds=time.monotonic()-started)
        atomic_json(out/'summary.json',summary)
        manifest={p.relative_to(out).as_posix():dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in sorted(out.rglob('*'))
                  if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json'] and 'null_runtime' not in p.parts}
        atomic_json(out/'output_manifest.json',manifest)
        state('complete',scientific_result=not args.smoke,files=len(manifest));atomic_json(out/'status.json',json.loads((runtime/'status.json').read_text()))
    except BaseException as error:
        state('failed',error=repr(error));(runtime/'traceback.txt').write_text(traceback.format_exc());raise
    finally:lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--smoke',action='store_true')
    args=parser.parse_args();settings=json.loads(args.config.read_text())
    with threadpool_limits(int(settings['threads'])):run(args)
