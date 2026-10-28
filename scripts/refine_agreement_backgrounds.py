"""Control-only refinement from saved profiles; preserve every scatter point.

Use the existing frozen, batch-matched public different-target/NTC plans.
Public NTC units are guide x batch, referenced to pooled NTC with whole-guide
exclusion, as in the common main methods. Compare only coverage-matched pairs.
"""
from pathlib import Path
import argparse,json,shutil
import numpy as np
import pandas as pd
from compare_agreement_overlay import members,atomic,sha,centers_from_means
from perturbseq_core_effects import pair_effects


def load_checkpoint(source,name):
    record=source/(name+'_records.csv');manifest=json.loads((source/'output_manifest.json').read_text())
    if sha(record)!=manifest[record.name]['sha256']:raise ValueError('Checkpoint record identity differs')
    records=pd.read_csv(record)
    # Older checkpoints wrote an object-typed guide label member. Never unpickle
    # it: the independently hashed companion CSV supplies the same saved order.
    profiles,proof=members(source/(name+'_profiles_checkpoint.npz'),['old','new','gene_ids'])
    if profiles['old'].shape!=profiles['new'].shape or profiles['old'].shape!=(len(records),len(profiles['gene_ids'])):raise ValueError('Checkpoint axes differ')
    profiles['guides']=records.guide.to_numpy(dtype=str)
    return profiles,records


def matched_coverage(plan,ratio=1.5):
    return ((plan.NTC_cells_a/plan.cells_a).between(1/ratio,ratio)
            &(plan.NTC_cells_b/plan.cells_b).between(1/ratio,ratio)
            &(plan.other_cells_a/plan.cells_a).between(1/ratio,ratio)
            &(plan.other_cells_b/plan.cells_b).between(1/ratio,ratio))


def pooled_control_profiles(controls,source,method):
    n=int(controls.cells.sum());weights=controls.cells.to_numpy()
    total=(source*weights[:,None]).sum(0);result=np.empty_like(source)
    own={g:np.flatnonzero(controls.guide.eq(g).to_numpy()) for g in controls.guide.unique()}
    for i,row in enumerate(controls.itertuples()):
        ix=own[row.guide];removed=weights[ix];remaining=n-int(removed.sum())
        if remaining<2:raise ValueError('Insufficient independent control reference')
        reference=(total-(source[ix]*removed[:,None]).sum(0))/remaining
        result[i]=np.log2((source[i]+1)/(reference+1)) if method=='old' else source[i]-reference
    return result


def public_control_effects(revision,genes):
    meta=pd.read_csv(revision/'tables/group_metadata.csv');nc=meta.source_control_type.eq('non_targeting_control').to_numpy()
    linear,lp=members(revision/'tables/A1_within_batch_mean_lfc_profiles.npz',['means','keys','gene_ids'])
    logged,bp=members(revision/'tables/B1_within_batch_cellwise_log2_profiles.npz',['means','keys','gene_ids'])
    if not np.array_equal(linear['keys'],logged['keys']) or not np.array_equal(meta.key,linear['keys']):raise ValueError('Public control row identities differ')
    if not np.array_equal(genes,linear['gene_ids']) or not np.array_equal(genes,logged['gene_ids']):raise ValueError('Public control RNA identities differ')
    controls=meta.loc[nc].reset_index(drop=True);n=int(controls.cells.sum())
    if n!=34632 or len(controls)!=624:raise ValueError('Public NTC cohort changed')
    effects={};arrays={}
    for method,source in [('old',linear['means'][nc]),('new',logged['means'][nc])]:
        effects[method]=pooled_control_profiles(controls,source,method)
    return controls,effects,dict(linear=lp,logged=bp)


def run(config):
    source=Path(config['source']);out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    revision=Path(config['public_revision']);baseline=Path(config['public_baseline'])
    summary=json.loads((source/'summary.json').read_text());common=pd.read_csv(source/'common_RNA.csv').gene_id.astype(str).to_numpy()
    paper_profiles,paper_records=load_checkpoint(source,'kolf');lookup={g:i for i,g in enumerate(paper_records.guide.astype(str))}
    ctrl,effects,control_proof=public_control_effects(revision,paper_profiles['gene_ids'])
    nc_lookup={g:i for i,g in enumerate(ctrl.key.astype(str))}
    frozen=pd.read_csv(baseline/'tables/diagnostic_matched_backgrounds.csv')
    original_different=frozen.loc[frozen.comparator.eq('different_target')].set_index(['guide_a','guide_b'])
    frozen_nc=pd.read_csv(revision/'tables/comparator_plan.csv')
    frozen_nc=frozen_nc.loc[frozen_nc.comparator.eq('within_batch_NTC_pair')&frozen_nc.matched.astype(bool)].set_index(['guide_a','guide_b'])
    plans={};counts={};bridge=[]
    for name in ['chip01','company','kolf']:
        plan=pd.read_csv(source/(name+'_comparison_plan.csv'))
        if name=='kolf':
            matched=[]
            for i,row in plan.iterrows():
                key=(row.guide_a,row.guide_b)
                d=original_different.loc[key]
                if not bool(d.matched):continue
                # Reproduce all old public different-target comparisons first.
                r=pair_effects(paper_profiles['old'][lookup[d.control_a]],paper_profiles['old'][lookup[d.control_b]])['pearson']
                bridge.append(abs(r-float(d.pearson)))
                if key not in frozen_nc.index:continue
                c=frozen_nc.loc[key]
                plan.loc[i,['other_a','other_b','NTC_a','NTC_b']]=[d.control_a,d.control_b,c.control_a,c.control_b]
                for prefix,a,b,rec in [('other',d.control_a,d.control_b,paper_records.set_index('guide')),
                                     ('NTC',c.control_a,c.control_b,ctrl.set_index('key'))]:
                    plan.loc[i,[prefix+'_cells_a',prefix+'_cells_b']]=[int(rec.loc[a,'cells']),int(rec.loc[b,'cells'])]
                matched.append(i)
            plan=plan.loc[matched]
        plan=plan.loc[matched_coverage(plan)].reset_index(drop=True)
        if not len(plan):raise ValueError('No coverage-matched controls for '+name)
        plans[name]=plan;counts[name]=dict(matched=len(plan),original=len(pd.read_csv(source/(name+'_comparison_plan.csv'))))
        plan.to_csv(out/(name+'_comparison_plan.csv'),index=False)
    if max(bridge)>1e-8:raise ValueError('Frozen public comparator bridge failed')
    for p in source.iterdir():
        if p.name.endswith('_points.csv') or p.name.endswith('_input_proof.json') or p.name=='common_RNA.csv':
            target=out/p.name;shutil.copyfile(p,target)
            if sha(p)!=sha(target):raise ValueError('Preserved point copy differs')
    for name,plan in plans.items():
        if name=='kolf':profiles=paper_profiles
        else:profiles,_=load_checkpoint(source,name)
        indices={g:i for i,g in enumerate(profiles['guides'].astype(str))}
        shared=pd.Index(profiles['gene_ids'].astype(str)).get_indexer(common)
        if (shared<0).any():raise ValueError('Common RNA lost')
        for suffix,selection in [('',slice(None)),('_commonRNA',shared)]:
            for method in ['old','new']:
                x=profiles[method][:,selection];nc=effects[method][:,selection] if name=='kolf' else x
                rows=[]
                for row in plan.itertuples():
                    comparisons=[('same_target',row.guide_a,row.guide_b,x,indices),('different_target',row.other_a,row.other_b,x,indices),('NTC',row.NTC_a,row.NTC_b,nc,nc_lookup if name=='kolf' else indices)]
                    for kind,a,b,array,index in comparisons:
                        result=pair_effects(array[index[a]],array[index[b]])
                        if not result['correlation_defined']:raise ValueError('Undefined fixed matched comparison')
                        rows.append(dict(method=method,target=row.target,comparison=kind,guide_a=a,guide_b=b,
                                         pearson=result['pearson'],features=result['features']))
                pd.DataFrame(rows).to_csv(out/f'{name}_{method}{suffix}_pairs.csv',index=False)
    summary.update(background_code=config['expected_commit'],background_matching=counts,
       public_old_different_target_bridge_max_error=max(bridge),
       public_NTC_unit='guide x batch; cell counts and batch composition matched; whole guide excluded from pooled reference',
       lab_NTC_unit='whole guide; nearest cell counts restricted to 1.5-fold; company capture labels unavailable',
       scatter_points_preserved=True,control_input_members=control_proof,
       inferential_pvalues=False,NTC_background='Fixed count-matched controls; known public batches matched. Descriptive, not FDR calibration.')
    atomic(out/'summary.json',summary)
    atomic(out/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()})
    print(json.dumps(dict(complete=True,matching=counts,old_public_comparator_error=max(bridge)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args();run(json.loads(args.config.read_text()))
