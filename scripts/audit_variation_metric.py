"""Frozen-cohort audit of the ten-cell average-absolute-distance comparison.

Uses saved guide means only. Full-cell L1, RMS distance, and response Pearson
are evaluated on the same RNA masks/pairs; no effect-dependent feature choice.
"""
from pathlib import Path
import argparse,io,json,os,time
import numpy as np
import pandas as pd
from compare_agreement_overlay import members,sha,atomic
from assess_screened_effects import reference_row


def metrics(a,b,ra,rb,keep):
    a,b,ra,rb=[np.asarray(x,float)[keep] for x in [a,b,ra,rb]]
    x,y=np.log2(a+1),np.log2(b+1)
    delta=x-y;ea=x-np.log2(ra+1);eb=y-np.log2(rb+1)
    r=float(np.corrcoef(ea,eb)[0,1]) if len(ea)>1 and ea.std()>1e-12 and eb.std()>1e-12 else None
    return dict(mean_absolute_difference=float(np.mean(abs(delta))),RMS_difference=float(np.sqrt(np.mean(delta**2))),response_pearson=r)


def sample_checks(plan,analysis_cells):
    parsed={}
    for row in plan.itertuples(index=False):
        values=np.fromstring(str(row.analysis_rows),sep=';',dtype=int)
        if len(values)!=int(row.cells) or len(set(values))!=len(values) or values.min()<0 or values.max()>=analysis_cells:
            raise ValueError('Invalid frozen sampled rows')
        key=(row.guide,int(row.cells),int(row.repeat),int(row.side))
        if key in parsed:raise ValueError('Duplicate sampling key')
        parsed[key]=set(values)
    for g,n,r,s in parsed:
        if s==0 and parsed[g,n,r,0]&parsed[g,n,r,1]:raise ValueError('Same-guide groups overlap')
        if n==5 and not parsed[g,5,r,s].issubset(parsed[g,10,r,s]):raise ValueError('5/10 cohorts not nested')
    return dict(operator_rows=len(parsed),disjoint_halves_checked=True,n5_nested_in_n10=True)


def run(config):
    out=Path(config['output']);out.mkdir(exist_ok=False)
    common=pd.read_csv(config['common_RNA']).gene_id.astype(str).to_numpy()
    if sha(Path(config['common_RNA']))!=config['common_RNA_sha256']:raise ValueError('Common RNA changed')
    started=time.time();summaries=[]
    for d in config['datasets']:
        dest=out/d['id'];dest.mkdir()
        atomic(out/'status.json',dict(stage='read_saved_profiles',dataset=d['id']))
        vr=Path(d['variation_root']);receipt=json.loads((vr/'complete.json').read_text())['files']
        for name in ['comparison_plan.csv','variation_points.csv','distance_draws.csv.gz','sampling_plan.csv.gz','eligible_guides.csv']:
            if sha(vr/name)!=receipt[name]['sha256']:raise ValueError('Accepted variation output changed')
        plan=pd.read_csv(vr/'comparison_plan.csv');old_points=pd.read_csv(vr/'variation_points.csv');draws=pd.read_csv(vr/'distance_draws.csv.gz')
        check=sample_checks(pd.read_csv(vr/'sampling_plan.csv.gz'),d['analysis_cells'])
        reproduced=draws.groupby(['cells_per_group','comparison','pair_id']).mean_absolute_difference.median()
        values=old_points.set_index(['cells_per_group','comparison','pair_id']).mean_absolute_difference
        if not np.allclose(reproduced.loc[values.index],values,atol=1e-12,rtol=1e-12):raise ValueError('Repeat aggregation differs')
        base=Path(d['baseline']);records=pd.read_csv(base/'tables/guide_effect_metadata.csv')
        if d.get('saved_profile_cache'):
            p=Path(d['saved_profile_cache']);r=json.loads(p.with_suffix('.json').read_text());body=p.read_bytes()
            import hashlib
            if hashlib.sha256(body).hexdigest()!=r['sha256']:raise ValueError('Accepted profile cache changed')
            with np.load(io.BytesIO(body),allow_pickle=False) as z:a={k:z[k] for k in z.files}
            source=dict(path=str(p),bytes=len(body),sha256=r['sha256'])
        else:
            a,source=members(base/'tables/rna_effect_profiles.npz',['means','gene_ids','guides'])
            i=int(np.flatnonzero(records.source_control_type.ne('non_targeting_control'))[0])
            a['center'],source['reference_prefix']=reference_row(base/'tables/rna_effect_profiles.npz',i)
            np.savez_compressed(dest/'saved_means.npz',**a)
            atomic(dest/'saved_means.json',dict(sha256=sha(dest/'saved_means.npz'),source=source))
        if not np.array_equal(a['guides'].astype(str),records.guide.astype(str)):raise ValueError('Guide profile axis differs')
        columns=pd.Index(a['gene_ids'].astype(str)).get_indexer(common)
        if (columns<0).any():raise ValueError('Common RNA absent')
        ratio=10000./d['baseline_scale'];means=a['means'][:,columns]*ratio;center=a['center'][columns]*ratio
        refs=np.broadcast_to(center,means.shape).copy()
        for i in np.flatnonzero(records.source_control_type.eq('non_targeting_control')):
            n=int(records.cells.iloc[i]);refs[i]=(center*d['NTC_cells']-n*means[i])/(d['NTC_cells']-n)
        lookup={g:i for i,g in enumerate(records.guide)};rna_index={g:i for i,g in enumerate(common)};direct={}
        for row in records.itertuples():
            if pd.notna(row.direct_feature_index):
                gene=str(a['gene_ids'][int(row.direct_feature_index)])
                if gene in rna_index:direct[row.guide]=rna_index[gene]
        all_rows=[];mask_sizes=[]
        for i,q in enumerate(plan.itertuples(index=False)):
            keep=np.ones(len(common),bool)
            for g in [q.guide_a,q.guide_b,q.other_a,q.other_b]:
                if g in direct:keep[direct[g]]=False
            mask_sizes.append(int(keep.sum()))
            for category,ga,gb in [('same_target',q.guide_a,q.guide_b),('different_target',q.other_a,q.other_b),('NTC',q.NTC_a,q.NTC_b)]:
                ia,ib=lookup[ga],lookup[gb]
                result=metrics(means[ia],means[ib],refs[ia],refs[ib],keep)
                all_rows.append(dict(dataset=d['id'],pair_id=i,target=q.target,comparison=category,guide_a=ga,guide_b=gb,
                                     full_cells_a=int(records.cells.iloc[ia]),full_cells_b=int(records.cells.iloc[ib]),RNA_features=int(keep.sum()),**result))
        masks=np.asarray(mask_sizes)
        for n in [5,10]:
            z=old_points.loc[old_points.cells_per_group.eq(n)]
            for category in z.comparison.unique():
                q=z.loc[z.comparison.eq(category)].sort_values('pair_id')
                if not np.array_equal(q.RNA_features.to_numpy(),masks):raise ValueError('Original RNA masks differ')
        full=pd.DataFrame(all_rows);full.to_csv(dest/'full_cell_pair_metrics.csv',index=False)
        subsets={'all':np.arange(len(plan))}
        if d.get('batch_patterns'):
            diagnostic=pd.read_csv(d['batch_patterns'])
            if not np.array_equal(diagnostic[['guide_a','guide_b','other_a','other_b']],plan[['guide_a','guide_b','other_a','other_b']]):raise ValueError('Batch pattern plan differs')
            subsets['exact_batch_pattern']=np.flatnonzero(diagnostic.exact_batch_pattern_match.eq(True))
        summary=[]
        for subset,ids in subsets.items():
            if not len(ids):continue
            for n in [5,10,'all_cells']:
                if n=='all_cells':table=full.loc[full.pair_id.isin(ids)];names=['mean_absolute_difference','RMS_difference','response_pearson']
                else:table=old_points.loc[old_points.pair_id.isin(ids)&old_points.cells_per_group.eq(n)];names=['mean_absolute_difference']
                for name in names:
                    z=table.pivot(index='pair_id',columns='comparison',values=name)
                    difference=z.different_target-z.same_target
                    advantage=-difference if name=='response_pearson' else difference
                    summary.append(dict(dataset=d['id'],subset=subset,cells_per_group=n,metric=name,pairs=len(z),
                                        same_target=float(z.same_target.median()),different_target=float(z.different_target.median()),
                                        NTC=float(z.NTC.median()),paired_advantage_median=float(advantage.median()),
                                        fraction_same_target_more_similar=float((advantage>0).mean())))
        canon=plan.apply(lambda r:'|'.join(sorted([r.NTC_a,r.NTC_b])),axis=1)
        ntc_dedup=full.loc[full.comparison.eq('NTC')].assign(NTC_pair=canon.to_numpy()).groupby('NTC_pair').mean_absolute_difference.median()
        summary_frame=pd.DataFrame(summary);summary_frame.to_csv(dest/'comparison_summary.csv',index=False)
        atomic(dest/'audit.json',dict(sample_checks=check,all_original_repeat_medians_reproduced=True,common_RNA=len(common),
                                     source=source,NTC_unique_pairs=int(canon.nunique()),most_reused_NTC_pair_occurrences=int(canon.value_counts().iloc[0]),
                                     full_NTC_unique_pair_weighted_median=float(ntc_dedup.median()),
                                     NTC_full_guide_batch_composition_not_matched=True,raw_count_sources_read=False,
                                     full_reference_profile_is_author_selected_strong=d['id']=='kolf'))
        summaries.extend(summary)
    pd.DataFrame(summaries).to_csv(out/'all_comparisons.csv',index=False)
    atomic(out/'status.json',dict(stage='complete',elapsed_seconds=time.time()-started,commit=config['expected_commit']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(json.loads(Path(p.parse_args().config).read_text()))
