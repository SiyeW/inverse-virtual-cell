"""Controlled repair of the public original-method challenge.

Reuse its exact QC cohort, RNA features, normalization and targeting draws.
First change only batch-matched references/control units; then compare the
established cell-wise log-before-mean representation on the same inputs.
Read only already accepted count-cache blocks; never the source H5AD payload.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import importlib.metadata as md
import io
import json
import os
from pathlib import Path
import signal
import time
import traceback

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu
from threadpoolctl import threadpool_limits

from perturbseq_core_effects import expression_effect, pair_effects, stable_rng, unique_direct_feature
from perturbseq_core_io import atomic_json, revision, sha256
from public_count_cache import axis_hash, check_arrays


def plan_groups(obs, baseline_records, seed=20261022, sample_sizes=(10, 50), repeats=8):
    """Targeting units stay identical; NTC guides are separated by batch."""
    obs = obs.copy()
    if not obs.index.is_unique or not obs.source_row.is_unique:
        raise ValueError('Nonunique current cohort identities')
    if not obs.n_gRNA.eq(1).all():
        raise ValueError('Expected adopted single-guide cohort')
    if obs.loc[obs.gene_target.ne('NTC')].groupby('guide_id').batch.nunique().max() != 1:
        raise ValueError('Targeting guide spans batches; cannot claim unchanged targeting unit')
    pos_by_guide = obs.groupby('guide_id', sort=True).indices
    ops = []; labels = []; starts = {}; columns = []; row_indices = []; weights = []
    def add(key, positions, **meta):
        positions = np.asarray(positions, dtype=np.int64)
        if len(positions) < 1:
            raise ValueError('Empty averaging unit')
        number = len(ops)
        row_indices.extend([number]*len(positions)); columns.extend(positions.tolist())
        weights.extend([1/len(positions)]*len(positions))
        ops.append(dict(key=key, cells=len(positions), **meta))
        return number
    starts['baseline'] = len(ops)
    for row in baseline_records.itertuples():
        positions = pos_by_guide[row.guide]
        if len(positions) != row.cells:
            raise ValueError('Original full-guide cell count changed')
        add(row.guide, positions)
    starts['groups'] = len(ops)
    groups = []
    grouped = obs.groupby(['batch', 'guide_id'], sort=True).indices
    for (batch, guide), positions in grouped.items():
        positions = np.asarray(positions)
        if len(positions) < 10:
            continue
        target = str(obs.gene_target.iloc[positions[0]])
        role = 'non_targeting_control' if target == 'NTC' else 'targeting'
        key = guide+'@'+batch if target == 'NTC' else guide
        row = dict(key=key, guide=guide, batch=batch, target=target, source_control_type=role,
                   cells=len(positions), positions=positions)
        row['mean_operator_row'] = add(key, positions)
        groups.append(row)
    for size in sample_sizes:
        starts['draws_'+str(size)] = len(ops)
        for group in groups:
            if group['cells'] < size:
                continue
            # The original guide key makes every targeting 10-cell draw identical.
            rng = stable_rng(seed, group['key'])
            for repeat in range(repeats):
                chosen = np.sort(rng.choice(group['positions'], size, replace=False))
                add(group['key'], chosen, kind='draw', n=size, repeat=repeat+1)
                labels.append(dict(key=group['key'], n=size, repeat=repeat+1,
                                   operator_row=len(ops)-1,
                                   source_rows=obs.iloc[chosen].source_row.astype(int).tolist()))
    starts['halves'] = len(ops)
    halves = []
    for group in groups:
        if group['cells'] < 20:
            continue
        positions = stable_rng(seed, 'halves:'+group['key']).permutation(group['positions'])
        for side, chosen in enumerate([positions[::2], positions[1::2]]):
            op = add(group['key'], chosen, kind='half', side=side)
            halves.append(dict(key=group['key'], side=side, operator_row=op, cells=len(chosen)))
    starts['references'] = len(ops)
    refs = []
    for batch in sorted(obs.batch.unique()):
        positions = np.flatnonzero((obs.batch.eq(batch) & obs.gene_target.eq('NTC')).to_numpy())
        if len(positions) < 100:
            raise ValueError('Too few within-batch NTC cells')
        refs.append(dict(batch=batch, operator_row=add(batch, positions), cells=len(positions)))
    operator = sparse.csr_matrix((weights, (row_indices, columns)), shape=(len(ops), len(obs)))
    np.testing.assert_allclose(np.asarray(operator.sum(1)).ravel(), 1., rtol=1e-12, atol=1e-12)
    meta = pd.DataFrame([{k:v for k,v in row.items() if k != 'positions'} for row in groups])
    return dict(operator=operator, groups=meta, draws=pd.DataFrame(labels),
                halves=pd.DataFrame(halves), references=pd.DataFrame(refs), starts=starts,
                sample_sizes=list(sample_sizes), repeats=repeats)


def stream_means(root, manifest, source_rows, gene_ids, obs, columns, operator, scale, progress):
    """One hash-checked byte read per consumed compressed block; bounded RAM."""
    root = Path(root)
    if (axis_hash(source_rows) != manifest['contract']['rows'] or
            axis_hash(gene_ids) != manifest['contract']['genes']):
        raise ValueError('Accepted source axes changed')
    columns = np.asarray(columns, dtype=np.int64)
    if (np.diff(columns) <= 0).any():
        raise ValueError('Fixed RNA columns must remain sorted unique')
    remap = np.full(len(source_rows), -1, dtype=np.int32)
    remap[obs.source_row.to_numpy(int)] = np.arange(len(obs), dtype=np.int32)
    factors = float(scale)/obs.total_counts.to_numpy(float)
    linear = np.empty((operator.shape[0], len(columns)), dtype=np.float64)
    logged = np.empty_like(linear)
    filled = np.zeros(len(columns), dtype=bool); consumed = 0
    required=[x for x in manifest['blocks'] if np.any((columns>=x['first'])&(columns<x['end']))]
    planned=sum(x['bytes'] for x in required)
    progress('cached_counts_to_fixed_profiles',completed_RNA=0,total_RNA=len(columns),
             cached_bytes_read=0,cached_bytes_planned=planned)
    for block in manifest['blocks']:
        use = np.flatnonzero((columns >= block['first']) & (columns < block['end']))
        if not len(use):
            continue
        path = root/block['file']
        before = path.stat()
        if before.st_size != block['bytes']:
            raise ValueError('Cached input size changed')
        payload = path.read_bytes()
        if hashlib.sha256(payload).hexdigest() != block['sha256']:
            raise ValueError('Cached input hash changed')
        with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
            data, rows, ptr = archive['data'], archive['rows'], archive['ptr']
        check_arrays(data, rows, ptr, len(source_rows))
        if (path.stat().st_size, path.stat().st_mtime_ns) != (before.st_size, before.st_mtime_ns):
            raise ValueError('Count cache changed while reading')
        vs = []; rr = []; pp = [0]
        for position in use:
            j = columns[position]-block['first']
            a,b = int(ptr[j]),int(ptr[j+1])
            r = remap[rows[a:b]]; keep = r >= 0
            kept = r[keep]
            values = data[a:b][keep].astype(float)*factors[kept]
            vs.append(values);rr.append(kept);pp.append(pp[-1]+len(kept))
        matrix = sparse.csc_matrix((np.concatenate(vs), np.concatenate(rr), pp),
                                    shape=(len(obs), len(use)))
        linear[:,use] = (operator @ matrix).toarray()
        matrix.data = np.log2(1+matrix.data)
        logged[:,use] = (operator @ matrix).toarray()
        filled[use] = True; consumed += len(payload)
        progress('cached_counts_to_fixed_profiles', completed_RNA=int(filled.sum()),
                 total_RNA=len(columns), cached_bytes_read=consumed,cached_bytes_planned=planned)
    if not filled.all():
        raise ValueError('Accepted cache does not cover every fixed original RNA')
    return linear, logged, consumed


def profiles_and_draws(averages, plan, selected_names, global_names, mode):
    groups = plan['groups'].copy()
    means = averages[groups.mean_operator_row.to_numpy(int)]
    refs = plan['references'].set_index('batch')
    centers = np.empty_like(means)
    for i,row in enumerate(groups.itertuples()):
        ref = refs.loc[row.batch]; center = averages[int(ref.operator_row)]
        count = int(ref.cells)
        if row.source_control_type == 'non_targeting_control':
            if count-row.cells < 2:
                raise ValueError('Within-batch leave-guide-out reference too small')
            center = (count*center-row.cells*means[i])/(count-row.cells)
            count -= row.cells
        centers[i] = center
        groups.loc[i,'reference_cells'] = count
    effects = np.log2((means+1)/(centers+1)) if mode == 'lfc' else means-centers
    direct = [unique_direct_feature(row.target, selected_names, global_names) for row in groups.itertuples()]
    for i,j in enumerate(direct):
        if j is not None: effects[i,j] = np.nan
    lookup = {key:i for i,key in enumerate(groups.key)}
    first = np.empty_like(effects); rows=[]
    for (key,n), draws in plan['draws'].groupby(['key','n'], sort=False):
        i=lookup[key]; values=averages[draws.operator_row.to_numpy(int)]
        changes=np.log2((values+1)/(centers[i]+1)) if mode=='lfc' else values-centers[i]
        if direct[i] is not None: changes[:,direct[i]]=np.nan
        for k,row in enumerate(draws.itertuples()):
            finite=changes[k,np.isfinite(changes[k])]
            rows.append(dict(key=key,guide=groups.guide.iloc[i],batch=groups.batch.iloc[i],
                source_control_type=groups.source_control_type.iloc[i],cells=int(n),repeat=int(row.repeat),
                rms=float(np.sqrt(np.mean(finite**2))),reference_cells=int(groups.reference_cells.iloc[i])))
            if n==10 and row.repeat==1:first[i]=changes[k]
    half_rows=[]
    for key,part in plan['halves'].groupby('key',sort=False):
        i=lookup[key];values=averages[part.sort_values('side').operator_row.to_numpy(int)]
        changes=np.log2((values+1)/(centers[i]+1)) if mode=='lfc' else values-centers[i]
        if direct[i] is not None:changes[:,direct[i]]=np.nan
        half_rows.append(dict(key=key,source_control_type=groups.source_control_type.iloc[i],
                              cells=int(groups.cells.iloc[i]),**pair_effects(changes[0],changes[1])))
    return groups,means,centers,effects,first,pd.DataFrame(rows),pd.DataFrame(half_rows)


def frozen_comparators(groups, original_pairs, original_background, seed=20261022):
    lookup=groups.set_index('key');rows=[]
    different=original_background.loc[original_background.comparator.eq('different_target')].set_index(['guide_a','guide_b'])
    ntc=groups.loc[groups.source_control_type.eq('non_targeting_control')]
    for pair in original_pairs.itertuples():
        a,b=lookup.loc[pair.guide_a],lookup.loc[pair.guide_b]
        chosen=different.loc[(pair.guide_a,pair.guide_b)]
        item=dict(guide_a=pair.guide_a,guide_b=pair.guide_b,target=pair.target,
                  comparator='different_target',matched=bool(chosen.matched))
        if item['matched']:item.update(control_a=chosen.control_a,control_b=chosen.control_b)
        rows.append(item)
        ca=ntc.loc[ntc.batch.eq(a.batch)&ntc.cells.between(a.cells/1.5,a.cells*1.5)]
        cb=ntc.loc[ntc.batch.eq(b.batch)&ntc.cells.between(b.cells/1.5,b.cells*1.5)]
        candidates=[(x.key,y.key) for x in ca.itertuples() for y in cb.itertuples() if x.guide!=y.guide]
        if candidates:
            rng=stable_rng(seed,'negative:'+pair.guide_a+':'+pair.guide_b)
            x,y=candidates[int(rng.integers(len(candidates)))]
            rows.append(dict(guide_a=pair.guide_a,guide_b=pair.guide_b,target=pair.target,
                             comparator='within_batch_NTC_pair',control_a=x,control_b=y,matched=True))
        else:
            rows.append(dict(guide_a=pair.guide_a,guide_b=pair.guide_b,target=pair.target,
                             comparator='within_batch_NTC_pair',matched=False))
    return pd.DataFrame(rows)


def evaluation(groups,effects,metrics,original_pairs,comparators):
    lookup={key:i for i,key in enumerate(groups.key)}
    pairs=[];background=[]
    for pair in original_pairs.itertuples():
        pairs.append(dict(target=pair.target,guide_a=pair.guide_a,guide_b=pair.guide_b,
                           **pair_effects(effects[lookup[pair.guide_a]],effects[lookup[pair.guide_b]])))
    for row in comparators.itertuples():
        if row.matched:
            background.append(dict(target=row.target,guide_a=row.guide_a,guide_b=row.guide_b,
                                   comparator=row.comparator,control_a=row.control_a,control_b=row.control_b,
                                   **pair_effects(effects[lookup[row.control_a]],effects[lookup[row.control_b]])))
    pairs=pd.DataFrame(pairs);background=pd.DataFrame(background)
    scores=metrics.groupby(['key','source_control_type','cells'],sort=True).rms.median().reset_index()
    scores=scores.merge(groups[['key','batch']],on='key',validate='many_to_one')
    comparison={}
    for n,part in scores.groupby('cells'):
        tc=part.loc[part.source_control_type.eq('targeting'),'rms'];ntc=part.loc[part.source_control_type.eq('non_targeting_control'),'rms']
        comparison[str(n)]=dict(targeting_guides=len(tc),NTC_guide_batch_units=len(ntc),
            targeting_RMS_median=float(tc.median()) if len(tc) else None,
            NTC_RMS_median=float(ntc.median()) if len(ntc) else None,
            targeting_vs_NTC_rank_fraction=float(mannwhitneyu(tc,ntc,method='asymptotic').statistic/(len(tc)*len(ntc))) if len(tc) and len(ntc) else None)
        comparison[str(n)]['by_batch']={}
        for batch,within in part.groupby('batch'):
            t=within.loc[within.source_control_type.eq('targeting'),'rms']
            c=within.loc[within.source_control_type.eq('non_targeting_control'),'rms']
            comparison[str(n)]['by_batch'][batch]=dict(targeting_units=len(t),NTC_units=len(c),
                targeting_RMS_median=float(t.median()) if len(t) else None,
                NTC_RMS_median=float(c.median()) if len(c) else None,
                targeting_vs_NTC_rank_fraction=float(mannwhitneyu(t,c,method='asymptotic').statistic/(len(t)*len(c))) if len(t) and len(c) else None)
    discrimination={}
    for kind,part in background.groupby('comparator'):
        joined=part.merge(pairs[['guide_a','guide_b','pearson']],on=['guide_a','guide_b'],suffixes=('_background','_same'),validate='one_to_one')
        difference=joined.pearson_same-joined.pearson_background
        discrimination[kind]=dict(matched_pairs=len(joined),all_original_pairs=len(pairs),
            unique_control_pairs=len({tuple(sorted([r.control_a,r.control_b])) for r in part.itertuples()}),
            same_pair_median=float(joined.pearson_same.median()),background_median=float(joined.pearson_background.median()),
            paired_r_advantage_median=float(difference.median()),fraction_same_above_background=float((difference>0).mean()))
    summary=dict(same_pair_median=float(pairs.pearson.median()),
        median_of_target_medians=float(pairs.groupby('target').pearson.median().median()),
        discrimination=discrimination,matched_RMS=comparison,
        inferential_pvalues='not_reported; guides, reused pairs and selected Strong are not independent validation units')
    # Compare sample sizes only on exactly the SAME guide x batch units.
    if 50 in set(scores.cells):
        wide=scores.pivot(index=['key','source_control_type'],columns='cells',values='rms').dropna(subset=[10,50])
        summary['precision_common50']={role:dict(units=len(part),RMS_n10_median=float(part[10].median()),
            RMS_n50_median=float(part[50].median())) for role,part in wide.groupby(level='source_control_type')}
    return summary,pairs,background


def revision_figures(out,source):
    """Derived diagnostic figures; no new filtering or scientific estimates."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    (out/'figures').mkdir()
    methods=['A0_original','A1_within_batch_mean_lfc','B1_within_batch_cellwise_log2']
    titles=['Original global reference','Within-batch reference','Within-batch, log before mean']
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),layout='constrained')
    for ax,name,title in zip(axes,methods,titles):
        if name=='A0_original':
            pair=pd.read_csv(source/'tables/within_target_pairs.csv')
            bg=pd.read_csv(source/'tables/diagnostic_matched_backgrounds.csv')
            choices=[('Different target',bg.loc[bg.comparator.eq('different_target')&bg.matched.astype(bool),'pearson'])]
        else:
            pair=pd.read_csv(out/f'tables/{name}_within_target_pairs.csv')
            bg=pd.read_csv(out/f'tables/{name}_background_pairs.csv')
            choices=[(kind,part.pearson) for kind,part in bg.groupby('comparator')]
        choices.insert(0,('Same target',pair.pearson))
        for label,values in choices:
            x=np.sort(pd.to_numeric(values,errors='raise').dropna().to_numpy())
            if len(x):ax.step(np.r_[x[0],x],np.r_[0,np.arange(1,len(x)+1)/len(x)],where='post',label=f'{label} ({len(x)})')
        ax.set(title=title,xlim=(-1,1),ylim=(0,1.02),xlabel='Guide-pair response Pearson r',ylabel='Cumulative fraction')
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Same original targeting pairs and RNA; matched backgrounds are required')
    fig.savefig(out/'figures/01_reference_specificity.png',dpi=160,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for ax,name,title in zip(axes,methods[1:],titles[1:]):
        table=pd.read_csv(out/f'tables/{name}_matched_RMS.csv')
        scores=table.groupby(['key','source_control_type','cells']).rms.median().reset_index()
        wide=scores.pivot(index=['key','source_control_type'],columns='cells',values='rms')
        if 50 not in wide.columns:ax.text(.5,.5,'No common 50-cell groups',ha='center');continue
        wide=wide.dropna(subset=[10,50])
        for role,color in [('non_targeting_control','#3c7c9f'),('targeting','#cc8940')]:
            part=wide.loc[wide.index.get_level_values('source_control_type')==role]
            for n,style in [(10,'--'),(50,'-')]:
                x=np.sort(part[n].to_numpy())
                if len(x):ax.step(np.r_[x[0],x],np.r_[0,np.arange(1,len(x)+1)/len(x)],where='post',
                    color=color,ls=style,label=f'{role}: n={n} ({len(x)})')
        ax.set(title=title,xlabel='RMS response (method-specific log2 units)',ylabel='Cumulative fraction',ylim=(0,1.02))
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Cell-number check on EXACTLY the same guide x batch units; 8 draws')
    fig.savefig(out/'figures/02_common_cohort_precision.png',dpi=160,bbox_inches='tight');plt.close(fig)
    cases=json.loads((out/'case_scatter_data.json').read_text())
    if cases:
        targets=sorted({x['target'] for x in cases})
        fig,axes=plt.subplots(len(targets),3,figsize=(15,4.5*len(targets)),squeeze=False,layout='constrained')
        for row,target in enumerate(targets):
            selected=[x for x in cases if x['target']==target]
            bound=max(max(abs(np.asarray(x[k]))) for x in selected for k in ['x','y'])*1.04
            for col,name in enumerate(methods):
                part=next(x for x in selected if x['method']==name);ax=axes[row,col]
                ax.scatter(part['x'],part['y'],s=3,alpha=.3)
                ax.plot([-bound,bound],[-bound,bound],ls='--',color='black',lw=.8)
                ax.set(xlim=(-bound,bound),ylim=(-bound,bound),xlabel=part['guide_a'],ylabel=part['guide_b'],
                       title=f'{target}: {titles[col]}\nr={part["pearson"]:.3f}; {len(part["x"])} RNA')
        fig.suptitle('Prespecified published reference cases; the original largest two guides')
        fig.savefig(out/'figures/03_published_reference_cases.png',dpi=160,bbox_inches='tight');plt.close(fig)


def run(config_path):
    config=json.loads(Path(config_path).read_text());code=revision(config['expected_commit'])
    source=Path(config['baseline']);out=Path(config['output'])
    if (source/'owner.lock').exists() or json.loads((source/'status.json').read_text())['stage']!='complete_pending_method_interpretation':
        raise ValueError('Completed immutable original baseline required')
    if not out.is_absolute() or out.parent!=Path(config['allowed_output_parent']):raise ValueError('New external result required')
    out.mkdir(exist_ok=False);(out/'tables').mkdir()
    owner=dict(pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],code=code)
    lock=out/'owner.lock';lock.write_text(json.dumps(owner));started=time.monotonic()
    def progress(stage,**more):
        status=dict(stage=stage,**owner,time=datetime.now().astimezone().isoformat(),elapsed_seconds=time.monotonic()-started,**more)
        atomic_json(out/'status.json',status);print(json.dumps(status),flush=True)
    def alarm(*_):raise TimeoutError('Protective revision deadline; original baseline and cache preserved')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(config['max_seconds'])
    try:
        with threadpool_limits(2):
            atomic_json(out/'config.resolved.json',config)
            progress('prepare_fixed_cohort_and_controls')
            published_manifest=json.loads((source/'output_manifest.json').read_text())
            original_manifest=published_manifest['outputs']
            if sha256(source/'raw_count_cache/manifest.json')!=published_manifest['raw_cache_manifest']:
                raise ValueError('Original count-cache manifest changed')
            for name in ['summary.json','matched_cell_selections.json','tables/all_released_cells_with_additional_QC.csv.gz',
                         'tables/all_RNA_feature_selection.csv.gz','tables/response_features.csv',
                         'tables/guide_effect_metadata.csv','tables/within_target_pairs.csv',
                         'tables/diagnostic_matched_backgrounds.csv']:
                if sha256(source/name)!=original_manifest[name]['sha256']:raise ValueError('Original metadata changed: '+name)
            baseline=json.loads((source/'summary.json').read_text())
            contract=baseline['contract']
            if any(contract[k]!=v for k,v in dict(raw_mean_gt=.5,matched_cells=10,matched_repeats=8,seed=20261022).items()):
                raise ValueError('The proposed controlled comparison does not match the original baseline')
            all_obs=pd.read_csv(source/'tables/all_released_cells_with_additional_QC.csv.gz',index_col=0)
            obs=all_obs.loc[all_obs.main_analysis_stage1.astype(bool)].copy()
            genes=pd.read_csv(source/'tables/all_RNA_feature_selection.csv.gz',index_col=0)
            feature=pd.read_csv(source/'tables/response_features.csv',index_col=0)
            columns=genes.index.get_indexer(feature.index)
            if (columns<0).any() or len(obs)!=baseline['additional_RNA_QC_cells']:raise ValueError('Original cohort/RNA changed')
            records=pd.read_csv(source/'tables/guide_effect_metadata.csv')
            pairs=pd.read_csv(source/'tables/within_target_pairs.csv')
            bg=pd.read_csv(source/'tables/diagnostic_matched_backgrounds.csv')
            plan=plan_groups(obs,records)
            # Enforce EXACT original targeting sample identities before reading cache.
            old_draws=json.loads((source/'matched_cell_selections.json').read_text())
            old_map={(x['guide'],x['repeat']):x['source_rows'] for x in old_draws}
            target_keys=set(plan['groups'].loc[plan['groups'].source_control_type.eq('targeting'),'key'])
            for row in plan['draws'].loc[plan['draws'].n.eq(10)&plan['draws'].key.isin(target_keys)].itertuples():
                if row.source_rows!=old_map[(row.key,row.repeat)]:raise ValueError('Targeting 10-cell draws changed')
            manifest=json.loads((source/'raw_count_cache/manifest.json').read_text())
            if manifest['contract']['source']!=baseline['source_identity']:raise ValueError('Cache input identity differs from accepted original')
            linear,logged,consumed=stream_means(source/'raw_count_cache',manifest,all_obs.index,genes.index,obs,columns,
                plan['operator'],baseline['normalization_target_umi'],progress)
            with np.load(source/'tables/rna_effect_profiles.npz',allow_pickle=False) as original:
                if not np.array_equal(original['gene_ids'].astype(str),feature.index.astype(str)):raise ValueError('Original feature axis differs')
                expected=original['means'];actual=linear[:len(records)]
                error=float(np.max(np.abs(expected-actual)))
                np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-10)
                first_target=int(np.flatnonzero(records.source_control_type.eq('targeting'))[0])
                old_global_reference=original['reference_means'][first_target].copy()
                original_effects=original['log2fc']
            offsets=[]
            for row in plan['references'].itertuples():
                center=linear[row.operator_row]
                delta=np.log2((old_global_reference+1)/(center+1))
                offsets.extend(dict(batch=row.batch,gene_id=str(gene),old_global_mean=float(old),
                                    within_batch_mean=float(new),additive_log2fc_offset=float(change))
                               for gene,old,new,change in zip(feature.index,old_global_reference,center,delta))
            pd.DataFrame(offsets).to_csv(out/'tables/reference_offsets.csv',index=False)
            atomic_json(out/'original_bridge.json',dict(all_guide_RNA_means_verified=True,max_absolute_error=error,
                original_targeting_draws_identical=True,cohort_cells=len(obs),RNA_features=len(feature),
                source_H5AD_payload_read=False,cached_bytes_read=consumed))
            plan['groups'].to_csv(out/'tables/group_metadata.csv',index=False)
            atomic_json(out/'sampling_plan.json',plan['draws'].to_dict('records'))
            comparators=frozen_comparators(plan['groups'],pairs,bg)
            comparators.to_csv(out/'tables/comparator_plan.csv',index=False)
            summaries={}
            plotted_effects={'A0_original':(records.guide.to_numpy(str),original_effects)}
            for name,averages,mode in [('A1_within_batch_mean_lfc',linear,'lfc'),('B1_within_batch_cellwise_log2',logged,'difference')]:
                progress('evaluate_revision',method=name)
                groups,means,centers,effects,first,metrics,halves=profiles_and_draws(averages,plan,feature.gene_name.to_numpy(str),genes.gene_name.to_numpy(str),mode)
                stats,pair_table,background=evaluation(groups,effects,metrics,pairs,comparators)
                stats['same_guide_half_medians']=halves.groupby('source_control_type').pearson.median().to_dict()
                stats['reference_cases']=json.loads(pair_table.loc[pair_table.target.isin(['RNF7','POU5F1'])].to_json(orient='records'))
                summaries[name]=stats
                plotted_effects[name]=(groups.key.to_numpy(str),effects)
                metrics.to_csv(out/f'tables/{name}_matched_RMS.csv',index=False)
                halves.to_csv(out/f'tables/{name}_same_guide_halves.csv',index=False)
                pair_table.to_csv(out/f'tables/{name}_within_target_pairs.csv',index=False)
                background.to_csv(out/f'tables/{name}_background_pairs.csv',index=False)
                np.savez_compressed(out/f'tables/{name}_profiles.npz',means=means,references=centers,effects=effects,
                    matched_first=first,keys=groups.key.to_numpy(str),gene_ids=feature.index.to_numpy(str))
            case_data=[]
            support=pd.read_csv(source/'tables/target_support.csv')
            for row in support.loc[support.target.isin(['RNF7','POU5F1'])].itertuples():
                for method,(keys,values) in plotted_effects.items():
                    lookup={g:i for i,g in enumerate(keys)}
                    if row.scatter_guide_a not in lookup or row.scatter_guide_b not in lookup:continue
                    x,y=values[lookup[row.scatter_guide_a]],values[lookup[row.scatter_guide_b]]
                    keep=np.isfinite(x)&np.isfinite(y)
                    case_data.append(dict(target=row.target,method=method,guide_a=row.scatter_guide_a,guide_b=row.scatter_guide_b,
                        pearson=pair_effects(x,y)['pearson'],gene_ids=feature.index.to_numpy(str)[keep].tolist(),
                        x=x[keep].tolist(),y=y[keep].tolist()))
            atomic_json(out/'case_scatter_data.json',case_data)
            summary=dict(code=code,baseline_code=baseline['code'],analysis_cells=len(obs),RNA_features=len(feature),
                source_count_payload_read=False,cached_bytes_read=consumed,original_bridge_max_error=error,
                changed_reference='within batch; whole own NTC guide excluded',
                matched_cells=[10,50],matched_repeats=8,seed=20261022,
                normalization_target_umi=baseline['normalization_target_umi'],
                targeting_unit_unchanged=True,NTC_unit='guide x batch; not independent biological replicates',
                methods=summaries,method_validation_complete=False,
                original_baseline_preserved=True,elapsed_seconds=time.monotonic()-started,
                versions={p:md.version(p) for p in ['numpy','pandas','scipy','h5py']})
            atomic_json(out/'summary.json',summary)
            progress('revision_diagnostic_figures')
            revision_figures(out,source)
        files={p.relative_to(out).as_posix():dict(bytes=p.stat().st_size,sha256=sha256(p))
            for p in out.rglob('*') if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json']}
        atomic_json(out/'output_manifest.json',files)
        progress('complete_pending_method_interpretation',files=len(files))
    except BaseException as error:
        progress('failed',error=repr(error));(out/'traceback.txt').write_text(traceback.format_exc());raise
    finally:
        signal.alarm(0);lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True)
    run(parser.parse_args().config)
