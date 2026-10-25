"""Finite NTC-reference sensitivity audit from accepted guide means, not cell counts."""
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import time
import zipfile

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from audit_rna_response_correlations import BudgetReader, corr, describe, identity, prefix_array, sha


def shift_stats(reference, original):
    delta = np.log2((reference + 1) / (original + 1))
    return dict(rms_log2_shift=float(np.sqrt(np.mean(delta**2))),
                median_abs_log2_shift=float(np.median(np.abs(delta))),
                p95_abs_log2_shift=float(np.quantile(np.abs(delta), .95)),
                max_abs_log2_shift=float(np.max(np.abs(delta))),
                genes_abs_shift_gt_005=int((np.abs(delta) > .05).sum()),
                genes_abs_shift_gt_01=int((np.abs(delta) > .1).sum()))


def calculate(records, means, original_reference, pairs, prefix_rows, splits=8, seed=20261025):
    ntc = records.index[records.source_control_type.eq('non_targeting_control')].to_numpy()
    if len(ntc) < 4 or int(ntc.max()) >= len(means) or len(means) != prefix_rows:
        raise ValueError('Incomplete whole-NTC guide coverage')
    counts = records.cells.to_numpy(float)
    n = counts[ntc]
    mu = means[ntc]
    reference = np.average(mu, axis=0, weights=n)
    np.testing.assert_allclose(reference, original_reference, atol=1e-10, rtol=1e-10)
    total = n.sum()
    rows, reference_rows, variants = [], [], []
    for j, i in enumerate(ntc):
        rest = (reference * total - mu[j] * n[j]) / (total - n[j])
        deviation = np.log2((mu[j] + 1) / (rest + 1))
        rows.append(dict(guide=records.iloc[i].guide, cells=int(n[j]), reference_cells=int(total-n[j]),
            reference_overlap_cells=0, cell_weight=float(n[j]/total),
            guide_rms_log2fc=float(np.sqrt(np.mean(deviation**2))),
            guide_median_log2fc=float(np.median(deviation)),
            **{'leave_one_'+k:v for k,v in shift_stats(rest, reference).items()}))
        reference_rows.append(dict(variant='leave_one_guide_out', name=str(records.iloc[i].guide),
                                   cells=int(total-n[j]), guides=len(ntc)-1, **shift_stats(rest,reference)))
    equal = mu.mean(axis=0)
    variants.append(('guide_equal_all',equal,equal,int(total),int(total),False))
    reference_rows.append(dict(variant='guide_equal_all',name='equal',cells=int(total),guides=len(ntc),**shift_stats(equal,reference)))
    split_rows=[]
    rng=np.random.default_rng(seed)
    for repeat in range(splits):
        order=rng.permutation(len(ntc));left,right=order[:len(order)//2],order[len(order)//2:]
        if set(left)&set(right) or set(left)|set(right)!=set(range(len(ntc))):
            raise ValueError('NTC partition overlap or missing groups')
        a,b=np.average(mu[left],axis=0,weights=n[left]),np.average(mu[right],axis=0,weights=n[right])
        na,nb=int(n[left].sum()),int(n[right].sum())
        label=str(repeat+1)
        for name,center,cells,members in [('a',a,na,left),('b',b,nb,right)]:
            reference_rows.append(dict(variant='whole_guide_half',name=label+name,cells=cells,guides=len(members),**shift_stats(center,reference)))
        split_rows.append(dict(repeat=repeat+1,cells_a=na,cells_b=nb,guides_a=len(left),guides_b=len(right),
            overlap_cells=0,ntc_guides_a=records.iloc[ntc[left]].guide.tolist(),ntc_guides_b=records.iloc[ntc[right]].guide.tolist(),
            **{'halves_'+k:v for k,v in shift_stats(a,b).items()}))
        variants.extend([(f'half_{label}_shared_a',a,a,na,na,False),
                         (f'half_{label}_shared_b',b,b,nb,nb,False),
                         (f'half_{label}_independent_ab',a,b,na,nb,True),
                         (f'half_{label}_independent_ba',b,a,nb,na,True)])
    lookup={g:i for i,g in enumerate(records.guide)}
    sensitivity=[]
    covered=pairs.loc[pairs.guide_a.map(lookup).lt(prefix_rows)&pairs.guide_b.map(lookup).lt(prefix_rows)]
    for pair in covered.itertuples():
        ia,ib=lookup[pair.guide_a],lookup[pair.guide_b]
        keep=np.ones(means.shape[1],bool)
        for i in [ia,ib]:
            value=records.iloc[i].direct_feature_index
            if pd.notna(value):keep[int(value)]=False
        base_a=np.log2((means[ia]+1)/(reference+1))
        base_b=np.log2((means[ib]+1)/(reference+1))
        original_r=corr(base_a[keep],base_b[keep])
        if original_r is None or abs(original_r-pair.pearson)>1e-10:
            raise ValueError('Accepted response pair differs')
        for label,ra,rb,na,nb,independent in [('cell_weighted_all',reference,reference,int(total),int(total),False)]+variants:
            a=np.log2((means[ia]+1)/(ra+1));b=np.log2((means[ib]+1)/(rb+1))
            value=corr(a[keep],b[keep])
            sensitivity.append(dict(target=pair.target,guide_a=pair.guide_a,guide_b=pair.guide_b,
                cells_a=pair.cells_a,cells_b=pair.cells_b,variant=label,independent_references=independent,
                reference_cells_a=na,reference_cells_b=nb,original_r=original_r,response_r=value,
                response_r_minus_original=None if value is None else value-original_r))
    guide_table=pd.DataFrame(rows);reference_table=pd.DataFrame(reference_rows);pair_table=pd.DataFrame(sensitivity)
    reference_by_gene=np.array([v[1] for v in variants if '_shared_' in v[0]])
    spread=np.log2((reference_by_gene+1)/(reference+1))
    gene_table=pd.DataFrame(dict(ntc_mean=reference,guide_equal_mean=equal,
        half_reference_log2_shift_sd=spread.std(axis=0,ddof=1),
        max_abs_half_reference_log2_shift=np.max(np.abs(spread),axis=0)))
    summary=dict(ntc_guides=len(ntc),ntc_cells=int(total),rna_features=means.shape[1],
        guide_cells=describe(n),largest_guide_cell_weight=float(n.max()/total),
        guide_weight_concentration_equivalent= float(1/np.sum((n/total)**2)),
        ntc_guide_rms=describe(guide_table.guide_rms_log2fc),
        max_leave_one_reference_rms=float(guide_table.leave_one_rms_log2_shift.max()),
        equal_guide_reference=shift_stats(equal,reference),
        half_reference_rms=describe(reference_table.loc[reference_table.variant.eq('whole_guide_half'),'rms_log2_shift']),
        split_reference_difference_rms=describe([x['halves_rms_log2_shift'] for x in split_rows]),
        reference_split_repeats=splits,seed=seed,covered_response_pairs=len(covered),covered_targets=int(covered.target.nunique()),
        full_response_pairs=len(pairs),pair_selection='all accepted same-target pairs fully contained in the required NTC row prefix; no correlation selection',
        common_reference_noise_note='For independent additive zero-mean errors, a shared noisy reference adds common covariance; it is not generally a mechanism forcing correlations downward.',
        individual_cell_variances_measured=False,true_biological_replicates=False,baseline_replacement_adopted=False,
        original_results_modified=False,new_p_values=False,per_cell_matrix_read=False)
    return guide_table,reference_table,pair_table,gene_table,split_rows,summary


def run(core,output,manifest_sha,commit,read_mib=8,seconds=360):
    started=time.monotonic();core=Path(core).resolve();output=Path(output).resolve()
    if output.exists() or output==core or output.is_relative_to(core) or core.is_relative_to(output):
        raise ValueError('Distinct unused audit directory required')
    if (core/'owner.lock').exists():raise ValueError('Active core writer')
    raw=(core/'output_manifest.json').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=manifest_sha:raise ValueError('Core manifest differs')
    manifest=json.loads(raw)
    for rel in ['tables/guide_effect_metadata.csv','tables/within_target_pairs.csv']:
        if sha(core/rel)!=manifest[rel]['sha256']:raise ValueError('Accepted small table differs')
    records=pd.read_csv(core/'tables/guide_effect_metadata.csv');pairs=pd.read_csv(core/'tables/within_target_pairs.csv')
    ntc=records.index[records.source_control_type.eq('non_targeting_control')].to_numpy()
    prefix_rows=int(ntc.max())+1
    if prefix_rows>256 or records.iloc[0].source_control_type=='non_targeting_control':
        raise ValueError('Unplanned prefix or invalid global-reference leading row')
    path=core/'tables/rna_effect_profiles.npz';before=identity(path)
    if before['bytes']!=manifest['tables/rna_effect_profiles.npz']['bytes']:raise ValueError('Profile size differs')
    print(json.dumps(dict(stage='reading_ntc_prefix',time=datetime.now().astimezone().isoformat(),prefix_rows=prefix_rows,
                         ntc_guides=len(ntc),read_budget_bytes=read_mib*1024**2,read_deadline_seconds=seconds)),flush=True)
    with path.open('rb') as stream:
        reader=BudgetReader(stream,read_mib*1024**2,time.monotonic()+seconds)
        with np.load(reader,allow_pickle=False) as values:guides,genes=values['guides'].astype(str),values['gene_ids'].astype(str)
        with zipfile.ZipFile(reader) as archive:
            means,shape=prefix_array(archive,'means',prefix_rows)
            reference,ref_shape=prefix_array(archive,'reference_means',1)
    if shape!=ref_shape or shape!=(len(records),len(genes)) or not np.array_equal(guides,records.guide.astype(str)):
        raise ValueError('RNA/guide axes differ')
    after=identity(path)
    if any(before[k]!=after[k] for k in ['bytes','mtime_ns','device']):raise ValueError('Source changed')
    print(json.dumps(dict(stage='input_complete',compressed_bytes=reader.bytes_read,elapsed_seconds=time.monotonic()-started)),flush=True)
    guide_table,reference_table,pair_table,gene_table,plans,summary=calculate(records,means,reference[0],pairs,prefix_rows)
    gene_table.insert(0,'gene_id',genes)
    summary.update(stage='complete',code_commit=commit,helper_sha256=sha(__file__),core_manifest_sha256=manifest_sha,
        source_identity_before=before,source_identity_after=after,source_npz_full_sha256_reused=manifest['tables/rna_effect_profiles.npz']['sha256'],
        source_fresh_full_hash=False,partial_numeric_member_crc_checked=False,global_reference_and_accepted_pair_parity_checked=True,
        compressed_profile_bytes_read=reader.bytes_read,read_budget_bytes=read_mib*1024**2,prefix_rows=prefix_rows,elapsed_seconds=time.monotonic()-started)
    output.mkdir()
    for name,table in [('ntc_guides',guide_table),('reference_variants',reference_table),('pair_sensitivity',pair_table),('gene_reference_variation',gene_table)]:
        table.to_csv(output/(name+'.csv'),index=False)
    (output/'reference_partitions.json').write_text(json.dumps(plans,indent=2),encoding='utf8')
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf8')
    proofs={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in output.iterdir()}
    (output/'output_manifest.json').write_text(json.dumps(proofs,indent=2),encoding='utf8')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['core','output','manifest-sha256','code-commit']:p.add_argument('--'+name,required=True)
    p.add_argument('--read-mib',type=int,default=8);p.add_argument('--read-seconds',type=int,default=360)
    a=p.parse_args()
    with threadpool_limits(limits=2):
        print(json.dumps(run(a.core,a.output,a.manifest_sha256,a.code_commit,a.read_mib,a.read_seconds),indent=2),flush=True)
