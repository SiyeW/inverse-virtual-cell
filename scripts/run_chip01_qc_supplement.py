"""Supplement fixed single-chip QC without changing counts or main analysis.

Private paths, identities, and bounded diagnostic budgets come from a local
configuration, never from defaults in this public script. No automatic rerun,
cell filtering, new clustering, or DNA-efficiency inference is performed.
"""
from pathlib import Path
from datetime import datetime
import argparse
import hashlib
import importlib.metadata as md
import json
import os
import subprocess
import time
import traceback

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

from analyze_expression_structure import sha256
from copy_verified_delivery import atomic_json
from run_chip01_literature import read_input, validate_core_policy_identity
from replogle2022_statistics import ADCalibration, ad_statistics
from perturbation_statistics import differential_expression_block
from chip01_qc_supplement import (distribution_summary, validate_guide_coverage,
    plan_null_splits, library_normalize, masked_response_profiles,
    guide_agreement, null_test_summary)


def save_figure(fig, path):
    fig.savefig(path, dpi=170, facecolor='white')
    import matplotlib.pyplot as plt
    plt.close(fig)


def ecdf(ax, x, label, color):
    x = np.sort(np.asarray(x, dtype=float))
    ax.step(x, np.arange(1, len(x)+1)/len(x), where='post', label=label, color=color)


def coverage_and_cell_qc(obs, coverage, output):
    import matplotlib.pyplot as plt
    table = validate_guide_coverage(obs, coverage)
    table.to_csv(output/'tables/guide_coverage.csv', index=False)
    summaries = []
    for phase, column in [('pre_RNA_QC', 'single_guide_before_qc'),
                          ('post_RNA_QC', 'single_guide_after_stage1_qc')]:
        for label, part in table.groupby('guide_type'):
            summaries.append(dict(phase=phase, guide_type=label, unit='library_guide',
                                  **distribution_summary(part[column], [1, 5, 10, 20, 50])))
    pd.DataFrame(summaries).to_csv(output/'tables/guide_coverage_summary.csv', index=False)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout='constrained')
    colors = {'NTC': '#3c7c9f', 'targeting': '#cc8940'}
    labels = ['NTC', 'targeting']
    for ax, (phase, column) in zip(axes[0], [('Before RNA QC', 'single_guide_before_qc'),
                                          ('After RNA QC', 'single_guide_after_stage1_qc')]):
        values = [table.loc[table.guide_type.eq(k), column].to_numpy() for k in labels]
        boxes = ax.boxplot(values, tick_labels=[f'{k}\n{len(v)} library guides' for k,v in zip(labels,values)],
                           patch_artist=True, showfliers=True, flierprops={'markersize': 2})
        for box, label in zip(boxes['boxes'], labels): box.set_facecolor(colors[label])
        ax.set(title=phase, ylabel='Single-guide cells per library guide', yscale='symlog')
    for label in labels:
        ecdf(axes[1,0], table.loc[table.guide_type.eq(label), 'single_guide_after_stage1_qc'], label, colors[label])
    axes[1,0].set(xscale='symlog', xlabel='Cells per guide after RNA QC (zeros retained)', ylabel='Fraction of library guides', title='Coverage ECDF')
    axes[1,0].legend(frameon=False)
    cutoffs = [1, 5, 10, 20]
    for offset, label in zip([-.18,.18], labels):
        x = table.loc[table.guide_type.eq(label), 'single_guide_after_stage1_qc'].to_numpy()
        axes[1,1].bar(np.arange(4)+offset, [100*np.mean(x < k) for k in cutoffs], width=.35, color=colors[label], label=label)
    axes[1,1].set(xticks=np.arange(4), xticklabels=['0 cells', '<5 cells', '<10 cells', '<20 cells'],
                  ylabel='Within-type library guides (%)', title='Zero/low coverage, not read depth')
    axes[1,1].legend(frameon=False)
    save_figure(fig, output/'figures/01_guide_coverage.png')

    single = obs.guide_status.astype(str).eq('single_guide')
    cell = obs.loc[single].copy()
    cell['guide_type'] = np.where(cell.is_control, 'NTC', 'targeting')
    metrics = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'raw_guide_umis', 'guide_umis']
    names = ['RNA UMIs', 'Detected RNA features', 'Mitochondrial RNA (%)',
             'All captured guide UMIs', 'Assigned guide UMIs']
    cell_summary = []
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    for ax, metric, name in zip(axes.flat, metrics, names):
        values, tick_labels = [], []
        for post in [False, True]:
            pool = cell.loc[cell.qc_pass_stage1] if post else cell
            for label in labels:
                x = pd.to_numeric(pool.loc[pool.guide_type.eq(label), metric]).to_numpy(dtype=float)
                values.append(x);tick_labels.append(label+('\npost' if post else '\npre'))
                cell_summary.append(dict(phase='post_RNA_QC' if post else 'pre_RNA_QC',
                                         guide_type=label, metric=metric, unit='single_guide_cell',
                                         **distribution_summary(x)))
        boxes = ax.boxplot(values, tick_labels=tick_labels, patch_artist=True,
                           showfliers=False, widths=.6)
        for box, label in zip(boxes['boxes'], labels*2): box.set_facecolor(colors[label])
        ax.set(title=name, ylabel=name)
        if metric != 'pct_counts_mt': ax.set_yscale('symlog')
    ax = axes[1,2]
    retained = []
    for label in labels:
        pool = cell.loc[cell.guide_type.eq(label)]
        retained.append(100*pool.qc_pass_stage1.mean())
    ax.bar(labels, retained, color=[colors[k] for k in labels])
    ax.set(ylabel='RNA-QC retained cells (%)', ylim=(0,100), title='Within single-guide calls')
    for i, value in enumerate(retained): ax.text(i, value, f'{value:.1f}%', ha='center', va='bottom')
    fig.suptitle('NTC versus targeting: pre/post RNA QC; boxplots omit fliers only for display')
    save_figure(fig, output/'figures/02_cell_quality.png')
    pd.DataFrame(cell_summary).to_csv(output/'tables/cell_quality_summary.csv', index=False)
    cell[[*metrics, 'guide_id', 'perturbed_gene', 'guide_type', 'qc_pass_stage1']].to_csv(
        output/'tables/single_guide_cell_qc.csv.gz')
    return table, pd.DataFrame(summaries), pd.DataFrame(cell_summary)


def response_agreement(counts, obs, genes, coverage, core, config, output):
    import matplotlib.pyplot as plt
    features = np.flatnonzero((genes.raw_mean.to_numpy() > .5) & genes.relative_z_defined.to_numpy(dtype=bool))
    if len(features) != config['expected_energy_features']:
        raise ValueError('Frozen response-feature identity differs')
    reference_total = float(np.median(obs.loc[core, 'total_counts']))
    norm = library_normalize(counts[:, features], obs.total_counts.to_numpy(), reference_total)
    records, profiles = masked_response_profiles(norm, obs.guide_id.astype(str), obs.perturbed_gene.astype(str),
        genes.core_control_normalized_mean.to_numpy()[features],
        genes.core_control_normalized_sd.to_numpy()[features],
        genes.gene_name.astype('string').fillna('').to_numpy(dtype=str)[features],
        genes.gene_name.astype('string').fillna('').to_numpy(dtype=str),
        min_cells=config['min_guide_cells'])
    return write_response_agreement(records, profiles, obs, coverage, genes, features,
                                    reference_total, config, output)


def write_response_agreement(records, profiles, obs, coverage, genes, features,
                             reference_total, config, output):
    import matplotlib.pyplot as plt
    pairs, valid = guide_agreement(records, profiles, max_between=config['max_between_pairs'], seed=config['seed'])
    records['profile_defined'] = valid
    records.to_csv(output/'tables/guide_response_profiles_metadata.csv', index=False)
    pairs.to_csv(output/'tables/guide_pair_correlations.csv.gz', index=False)
    np.savez_compressed(output/'tables/guide_response_profiles.npz', values=profiles,
        guides=records.guide.to_numpy(dtype=str), targets=records.target.to_numpy(dtype=str),
        feature_indices=features, gene_ids=genes.index.to_numpy(dtype=str)[features])
    targets = []
    for target in sorted(set(coverage.loc[coverage.guide_type.eq('targeting'), 'target_gene_name'].astype(str))):
        found = records.loc[records.target.eq(target) & records.profile_defined]
        source = obs.loc[obs.perturbed_gene.astype(str).eq(target)]
        within = pairs.loc[pairs.pair_type.eq('same-target') & pairs.target_a.eq(target), 'response_correlation']
        targets.append(dict(target=target, cells=len(source), observed_guides=source.guide_id.nunique(),
            eligible_guides=len(found), eligible_cells=int(found.cells.sum()), guide_pairs=len(within),
            median_correlation=float(within.median()) if len(within) else None,
            min_correlation=float(within.min()) if len(within) else None,
            max_correlation=float(within.max()) if len(within) else None,
            assessment='descriptive_assessable' if len(within) else 'insufficient_two_guide_profiles'))
    targets = pd.DataFrame(targets)
    targets.to_csv(output/'tables/target_guide_agreement.csv', index=False)
    kinds = ['same-target', 'different-target', 'NTC-NTC']
    colors = ['#3c7c9f', '#c89255', '#8b9298']
    fig, axes = plt.subplots(1,2,figsize=(12,5),layout='constrained')
    values = [pairs.loc[pairs.pair_type.eq(k), 'response_correlation'].to_numpy() for k in kinds]
    boxes = axes[0].boxplot(values, tick_labels=[f'{k}\nn={len(v)} pairs' for k,v in zip(kinds,values)],
                            patch_artist=True, showfliers=False)
    for box, color in zip(boxes['boxes'], colors):box.set_facecolor(color)
    axes[0].axhline(0,color='grey',lw=.8)
    axes[0].set(ylabel='Mean-z RNA response Pearson r', title=f'Each guide >= {config["min_guide_cells"]} cells')
    for kind,color in zip(kinds,colors):
        part=pairs.loc[pairs.pair_type.eq(kind)]
        if len(part)>1500:part=part.sample(1500,random_state=config['seed'])
        axes[1].scatter(part.min_cells,part.response_correlation,s=9,alpha=.4,color=color,label=kind)
    axes[1].set(xlabel='Smaller cell count in guide pair',ylabel='Response correlation',xscale='log', title='Coverage remains a limitation')
    axes[1].legend(frameon=False,fontsize=8)
    fig.suptitle('Descriptive guide agreement; guide pairs are not biological replicates')
    save_figure(fig,output/'figures/03_guide_agreement.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    bins=np.minimum(targets.eligible_guides.to_numpy(),3)
    axes[0].bar(np.arange(4),[(bins==i).sum() for i in range(4)],color='#6f99ad')
    axes[0].set(xticks=np.arange(4),xticklabels=['0','1','2','>=3'],xlabel='Eligible guide profiles per targeting entity',ylabel='Targeting entities')
    part=targets.loc[targets.guide_pairs>0]
    axes[1].scatter(part.cells,part.median_correlation,s=25,color='#3c7c9f',alpha=.7)
    axes[1].axhline(0,color='grey',lw=.8)
    axes[1].set(xlabel='Total RNA-QC single-guide cells per target',ylabel='Median within-target response r',title='Only assessable targets')
    save_figure(fig,output/'figures/04_target_agreement_coverage.png')
    summary = dict(feature_count=len(features),min_guide_cells=config['min_guide_cells'],
        profile_guides=len(records),profile_undefined=int((~valid).sum()),targeting_library_targets=len(targets),
        targets_assessable=int((targets.guide_pairs>0).sum()),
        targets_unassessable=int((targets.guide_pairs==0).sum()),
        pair_counts={k:len(v) for k,v in zip(kinds,values)},
        pair_median_r={k:float(np.median(v)) if len(v) else None for k,v in zip(kinds,values)},
        between_pair_subsample_cap=config['max_between_pairs'],
        normalization_target=reference_total,
        interpretation='Descriptive; shared baseline, sample sizes, pair dependence and single-chip design prevent biological replication claims')
    atomic_json(output/'guide_agreement_summary.json',summary)
    return summary


def negative_control_checks(counts, obs, genes, core, config, output, runtime, state, smoke):
    import matplotlib.pyplot as plt
    indices=np.flatnonzero(genes.relative_z_defined.to_numpy(dtype=bool))
    if len(indices)!=config['expected_tested_genes']:raise ValueError('Frozen DE testing domain differs')
    ntc=np.flatnonzero(obs.is_control.to_numpy(dtype=bool))
    ref_total=float(config.get('normalization_target_umi',np.median(obs.loc[core,'total_counts'])))
    norm=library_normalize(counts[ntc],obs.iloc[ntc].total_counts.to_numpy(),ref_total)
    sizes=config['null_sizes'] if not smoke else [1,25]
    repeats=config['null_repeats'] if not smoke else 2
    n_ref=config['null_reference_cells'] if not smoke else 200
    if smoke:indices=indices[np.linspace(0,len(indices)-1,min(128,len(indices)),dtype=int)]
    splits=plan_null_splits(len(ntc),sizes,n_ref,repeats,config['seed'])
    selection={s['comparison']:dict(treatment_barcodes=obs.index[ntc[s['treatment']]].astype(str).tolist(),
                 control_barcodes=obs.index[ntc[s['control']]].astype(str).tolist(),cells=s['cells']) for s in splits}
    atomic_json(output/'null_selection.json',selection)
    ad_p=np.lib.format.open_memmap(runtime/'ad_pvalues.npy',mode='w+',dtype=float,shape=(len(splits),len(indices)))
    mw_p=np.lib.format.open_memmap(runtime/'mw_pvalues.npy',mode='w+',dtype=float,shape=ad_p.shape)
    ad_p[:]=np.nan;mw_p[:]=np.nan
    detection=np.empty(len(indices));calibration=ADCalibration(Path(config['calibration']))
    started=time.monotonic();block_size=config['null_gene_block_size']
    total_blocks=(len(indices)+block_size-1)//block_size
    for block,start in enumerate(range(0,len(indices),block_size)):
        end=min(start+block_size,len(indices))
        dense=norm[:,indices[start:end]].toarray()
        detection[start:end]=(dense>0).mean(0)
        for row,split in enumerate(splits):
            x=dense[split['treatment']];control=dense[split['control']]
            xs=np.asfortranarray(np.sort(x,axis=0));ys=np.asfortranarray(np.sort(control,axis=0))
            ad_p[row,start:end]=calibration.pvalues(ad_statistics(xs,ys,already_sorted=True))
            if split['cells']>=10:
                mw_p[row,start:end]=differential_expression_block(x,control)['pvalue']
        ad_p.flush();mw_p.flush()
        state('null_gene_blocks',completed_blocks=block+1,total_blocks=total_blocks,
              completed_genes=end,total_genes=len(indices),comparisons=len(splits),
              stage_elapsed_seconds=time.monotonic()-started,
              seconds_per_block=(time.monotonic()-started)/(block+1))
    rows=[]
    for i,split in enumerate(splits):
        for method,values in [('AD',ad_p[i]),('MW',mw_p[i])]:
            if method=='MW' and split['cells']<10:continue
            rows.append(dict(comparison=split['comparison'],cells=split['cells'],replicate=split['replicate'],
                             method=method,control_cells=n_ref,**null_test_summary(values)))
    results=pd.DataFrame(rows);results.to_csv(output/'tables/null_comparison_summary.csv',index=False)
    np.savez_compressed(output/'tables/null_pvalues.npz',ad_pvalues=ad_p,mw_pvalues=mw_p,
        comparisons=np.array([s['comparison'] for s in splits]),sizes=np.array([s['cells'] for s in splits]),
        gene_ids=genes.index.to_numpy(dtype=str)[indices],ntc_fraction_detected=detection)
    summary=results.groupby(['method','cells']).agg(comparisons=('comparison','size'),
        median_raw_p_lt_005=('raw_p_lt_005','median'),median_bh_genes=('bh_q_lt_005','median'),
        maximum_bh_genes=('bh_q_lt_005','max'),contrasts_with_any_bh=('bh_q_lt_005',lambda x:int((x>0).sum()))).reset_index()
    summary['fraction_contrasts_with_any_bh']=summary.contrasts_with_any_bh/summary.comparisons
    summary.to_csv(output/'tables/null_by_sample_size.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    palette=plt.cm.viridis(np.linspace(.05,.9,len(sizes)))
    for ax,method,values in zip(axes,['AD','MW'],[ad_p,mw_p]):
        xmax=0
        for n,color in zip(sizes,palette):
            if method=='MW' and n<10:continue
            use=np.array([s['cells']==n for s in splits])
            p=np.sort(np.asarray(values[use]).ravel());p=p[np.isfinite(p)]
            ranks=np.unique(np.geomspace(1,len(p),300).astype(int))-1
            expected=-np.log10((ranks+.5)/len(p));observed=-np.log10(np.maximum(p[ranks],1e-300))
            xmax=max(xmax,float(expected.max()))
            ax.plot(expected,observed,label=f'n={n}',color=color,lw=1.2)
        ax.plot([0,xmax],[0,xmax],ls='--',color='black',lw=1,label='Uniform reference')
        ax.set(xlabel='Expected -log10(p), uniform reference',ylabel='Observed -log10(p)',title=method+' negative-control p-values')
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Randomly labelled original NTCs; discrete tests may be conservative')
    save_figure(fig,output/'figures/05_null_pvalue_qq.png')
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,method,color in zip(axes,['AD','MW'],['#bd7640','#3c7c9f']):
        part=results.loc[results.method.eq(method)]
        for n in sorted(part.cells.unique()):
            values=part.loc[part.cells.eq(n),'bh_q_lt_005'].to_numpy()
            offsets=np.linspace(-.08,.08,len(values))
            ax.scatter(np.full(len(values),n)*(1+offsets),values,s=30,color=color,alpha=.65)
        ax.set(xscale='log',yscale='symlog',xlabel='Pseudo-perturbation cells (NTC only)',
               ylabel='Genes with full-family BH q < 0.05',title=method+' calls on randomized negative controls')
        ax.axhline(0,color='grey',lw=.8)
    fig.suptitle('Each point is one pseudo-contrast; not independent biological replication')
    save_figure(fig,output/'figures/06_null_bh_calls.png')
    out=dict(original_ntc_cells=len(ntc),reference_cells_per_contrast=n_ref,
        pseudo_comparisons=len(splits),tested_genes=len(indices),sizes=sizes,repeats_per_size=repeats,
        seed=config['seed'],internal_overlap=False,full_family_bh=True,
        estimates=summary.to_dict(orient='records'),
        interpretation='Within-chip diagnostic on random NTC labels; no recalibration of original DE, no guarantee of biological FDR; repeated contrasts reuse cells')
    atomic_json(output/'negative_control_summary.json',out)
    return out


def run(args):
    config=json.loads(args.config.read_text())
    output=Path(config['output']);runtime=Path(config['runtime'])
    source=Path(config['source'])
    if output.parent.resolve()!=Path(config['allowed_output_parent']).resolve():
        raise ValueError('Output outside exact approved parent')
    if output.exists() or runtime.exists():raise FileExistsError('Distinct output and runtime required; inspect before any recovery')
    if output.resolve() in source.resolve().parents or source.resolve() in output.resolve().parents:
        raise ValueError('Source/output overlap')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip()
    if commit!=config['expected_commit']:raise ValueError('Execution code is not the pinned pushed commit')
    output.mkdir();runtime.mkdir(parents=True)
    for name in ['tables','figures']: (output/name).mkdir()
    owner=dict(pid=os.getpid(),uid=os.getuid(),code_commit=commit,started_at=datetime.now().astimezone().isoformat())
    lock=output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    with os.fdopen(fd,'w') as f:json.dump(owner,f)
    started=time.monotonic()
    def state(stage,**extra):
        value=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**owner,
                   elapsed_seconds=time.monotonic()-started)
        value.update(extra)
        atomic_json(runtime/'status.json',value)
        print(json.dumps(value),flush=True)
    source_before=source.stat()
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
        state('verify_selected_inputs')
        inputs={}
        for key in ['coverage','core_policy','core_guides','gene_statistics','calibration']:
            p=Path(config[key]);inputs[key]=dict(path=str(p),bytes=p.stat().st_size,sha256=sha256(p))
        if inputs['calibration']['sha256']!=config['calibration_sha256']:raise ValueError('Frozen AD reference differs')
        atomic_json(output/'input_manifest.json',inputs)
        state('verify_source_hash_and_read_sparse_counts')
        counts,obs,var,identity=read_input(source,config['source_sha256'])
        if identity['source_shape']!=config['expected_shape']:raise ValueError('Source dimensions differ')
        policy=pd.read_csv(config['core_policy'],index_col=0)
        validate_core_policy_identity(policy,obs)
        core=policy.core_control.to_numpy(dtype=bool)
        if int(core.sum())!=config['expected_core_cells']:raise ValueError('Core-reference cell identity differs')
        genes=pd.read_csv(config['gene_statistics'],index_col=0)
        if not np.array_equal(genes.index.astype(str).to_numpy(),var.index.astype(str).to_numpy()):raise ValueError('Gene row identity/order differs')
        source_view=ad.read_h5ad(source,backed='r');all_obs=source_view.obs.copy();source_view.file.close()
        coverage=pd.read_csv(config['coverage'])
        state('stratified_qc')
        reference,guide_summary,cell_summary=coverage_and_cell_qc(all_obs,coverage,output)
        state('guide_response_agreement')
        agreement=response_agreement(counts,obs,genes,reference,core,config,output)
        state('negative_control_checks')
        null=negative_control_checks(counts,obs,genes,core,config,output,runtime,state,args.smoke)
        state('verify_outputs_and_preserved_input_identity')
        after=source.stat()
        if (source_before.st_size,source_before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise ValueError('Source identity changed during read-only analysis')
        for key,entry in inputs.items():
            if sha256(Path(entry['path']))!=entry['sha256']:raise ValueError('Selected source artifact changed: '+key)
        summary=dict(run_id=config['run_id'],code_commit=commit,scientific_result=not args.smoke,
            engineering_smoke=bool(args.smoke),source_identity=identity,source_modified=False,
            library_guides=len(reference),ntc_library_guides=int(reference.guide_type.eq('NTC').sum()),
            observed_ntc_guides_after_qc=int(((reference.guide_type=='NTC')&(reference.single_guide_after_stage1_qc>0)).sum()),
            original_ntc_cells=int(obs.is_control.sum()),core_ntc_cells=int(core.sum()),
            guide_agreement=agreement,negative_controls=null,
            cell_filtering_changed=False,new_clustering_run=False,dna_editing_efficiency='not_assessed',
            versions={name:md.version(name) for name in ['numpy','scipy','pandas','anndata','numba','matplotlib']},
            settings=config,elapsed_seconds=time.monotonic()-started)
        atomic_json(output/'summary.json',summary)
        manifest={}
        for p in sorted(output.rglob('*')):
            if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json']:
                manifest[p.relative_to(output).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha256(p))
        atomic_json(output/'output_manifest.json',manifest)
        state('complete',scientific_result=not args.smoke,files=len(manifest))
        atomic_json(output/'status.json',json.loads((runtime/'status.json').read_text()))
    except BaseException as error:
        state('failed',error=repr(error))
        (runtime/'traceback.txt').write_text(traceback.format_exc())
        raise
    finally:
        lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--smoke',action='store_true')
    args=parser.parse_args()
    settings=json.loads(args.config.read_text())
    with threadpool_limits(int(settings['threads'])):
        run(args)
