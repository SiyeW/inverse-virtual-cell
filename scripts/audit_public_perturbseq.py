"""Bounded count-based audit of an already processed public Perturb-seq H5AD.

No dense X loading, new cell filtering, hit calling, or biological replicate
inference. Equal-sized disjoint arms and whole-guide held-out NTC references
make descriptive correlations comparable. Selection of the distributed cohort
remains an external ascertainment limitation.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import h5py
import numpy as np
import pandas as pd


def read_field(node):
    if isinstance(node, h5py.Group):
        codes = node['codes'][:]
        if (codes < 0).any():
            raise ValueError('Required metadata has missing categorical codes')
        return node['categories'].asstr()[:][codes]
    return node.asstr()[:] if node.dtype.kind in 'OSU' else node[:]


def stable_seed(seed, *labels):
    text = json.dumps([seed, *map(str, labels)], ensure_ascii=True)
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], 'little')


def make_plan(obs, arm_sizes, seed):
    """References are whole-guide disjoint from every test/control arm."""
    obs = obs.copy()
    if obs.groupby('gRNA').gene_target.nunique().max() != 1:
        raise ValueError('A guide maps to multiple targets')
    ntc = obs.gene_target.eq('NTC').to_numpy()
    guides = sorted(obs.loc[ntc, 'gRNA'].unique())
    if len(guides) < 12:
        raise ValueError('Too few independent NTC guide identities')
    rng = np.random.default_rng(seed)
    shuffled = np.asarray(guides)[rng.permutation(len(guides))]
    n_ref = max(2, len(guides) // 6)
    ref_a, ref_b = set(shuffled[:n_ref]), set(shuffled[n_ref:2*n_ref])
    assert ref_a.isdisjoint(ref_b)
    ref = np.full(len(obs), -1, dtype=np.int32)
    for bi, batch in enumerate(sorted(obs.batch.unique())):
        for half, names in enumerate([ref_a, ref_b]):
            rows = obs.batch.eq(batch) & obs.gRNA.isin(names)
            if rows.sum() < 100:
                raise ValueError('A batch lacks sufficient disjoint reference cells')
            ref[rows.to_numpy()] = bi * 2 + half
    records, maps = [], {}
    full = np.empty(len(obs), dtype=np.int32)
    for (batch, guide), positions in obs.groupby(['batch', 'gRNA'], sort=True).indices.items():
        positions = np.asarray(positions)
        target = str(obs.gene_target.iloc[positions[0]])
        gid = len(records)
        full[positions] = gid
        records.append(dict(group=gid, batch=str(batch), guide=str(guide),
                            target=target, role='ntc' if target == 'NTC' else 'targeting',
                            cells=len(positions), reference_guide=guide in ref_a | ref_b))
    for n in arm_sizes:
        assignments = np.full(len(obs), -1, dtype=np.int32)
        entries = []
        for row in records:
            if row['reference_guide'] or row['cells'] < 2*n:
                continue
            positions = np.flatnonzero(full == row['group'])
            rng = np.random.default_rng(stable_seed(seed, row['batch'], row['guide']))
            positions = rng.permutation(positions)[:2*n]
            first = len(entries)
            assignments[positions[:n]], assignments[positions[n:]] = first, first+1
            for half in [0, 1]:
                entries.append(dict(**row, arm=half, arm_id=first+half, arm_cells=n))
        maps[int(n)] = dict(assignment=assignments, entries=pd.DataFrame(entries))
    return pd.DataFrame(records), full, ref, maps, sorted(ref_a), sorted(ref_b)


def read_counts(node, progress, chunk_entries=4*1024**2):
    """One sequential read of each sparse value/index dataset, never dense X."""
    encoding = node.attrs.get('encoding-type')
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding != 'csc_matrix':
        raise ValueError('Declared counts layer must be CSC')
    shape = tuple(map(int, node.attrs['shape']))
    ptr = node['indptr'][:].astype(np.int64)
    if len(ptr) != shape[1]+1 or ptr[0] != 0 or (np.diff(ptr)<0).any():
        raise ValueError('Invalid sparse column pointer')
    count = int(ptr[-1])
    if len(node['data']) != count or len(node['indices']) != count:
        raise ValueError('Sparse value/index cardinality differs')
    if count * 8 > 6 * 1024**3:
        raise ValueError('Sparse payload exceeds the declared 6 GiB read budget')
    data = np.empty(count, dtype=np.float32)
    indices = np.empty(count, dtype=np.int32)
    digests = {}
    total = data.nbytes + indices.nbytes
    complete = 0
    for name, target in [('data', data), ('indices', indices)]:
        digest = hashlib.sha256()
        for first in range(0, count, chunk_entries):
            end = min(first+chunk_entries, count)
            node[name].read_direct(target, np.s_[first:end], np.s_[first:end])
            block = target[first:end]
            if name == 'data':
                if not np.isfinite(block).all() or (block<=0).any() or not np.equal(block,np.floor(block)).all():
                    raise ValueError('Expected positive raw integer nonzero counts')
            elif (block<0).any() or (block>=shape[0]).any():
                raise ValueError('Sparse row identity out of bounds')
            digest.update(memoryview(block).cast('B'))
            complete += block.nbytes
            if first == 0 or end == count or complete // (256*1024**2) != (complete-block.nbytes)//(256*1024**2):
                progress('read_counts', logical_bytes=complete, logical_total=total)
        digests[name] = digest.hexdigest()
    totals = np.zeros(shape[0], dtype=np.float64)
    detected = np.zeros(shape[0], dtype=np.int64)
    for first in range(0, count, chunk_entries):
        end = min(first+chunk_entries, count)
        totals += np.bincount(indices[first:end],weights=data[first:end],minlength=shape[0])
        detected += np.bincount(indices[first:end],minlength=shape[0])
    for lo,hi in zip(ptr[:-1],ptr[1:]):
        if np.any(np.diff(indices[lo:hi])<=0):
            raise ValueError('Raw count columns must be sorted and duplicate-free')
    if (totals<=0).any():
        raise ValueError('Retained record has no RNA counts')
    return data, indices, ptr, totals, detected, digests


def aggregate_column(values, rows, assignment, n_groups):
    keys = assignment[rows]
    keep = keys >= 0
    return np.bincount(keys[keep],weights=values[keep],minlength=n_groups)


def select_panel(h, obs, panel_size=512, target_size=64):
    """Pre-existing HVG ranks and coverage, never observed effect/correlation."""
    names=np.asarray(read_field(h['var/_index']),dtype=str)
    hvg=h['var/highly_variable'][:].astype(bool)
    rank=h['var/highly_variable_rank'][:].astype(float)
    eligible=np.flatnonzero(hvg & np.isfinite(rank))
    ordered=eligible[np.lexsort((eligible,rank[eligible]))]
    if len(ordered)<panel_size:
        raise ValueError('Insufficient pre-existing HVG ranks for declared panel')
    rna=ordered[:panel_size]
    coverage=obs.loc[~obs.gene_target.eq('NTC')].groupby('gene_target',sort=True).size()
    # Stable tie handling: alphabetical target identity after descending coverage.
    ranked=coverage.rename('n').reset_index().sort_values(['n','gene_target'],ascending=[False,True])
    targets=ranked.loc[ranked.gene_target.isin(names),'gene_target'].head(target_size).tolist()
    extra=np.flatnonzero(np.isin(names,targets)|np.char.startswith(np.char.upper(names),'MT-'))
    selected=np.unique(np.concatenate([rna,extra]))
    ptr=h['layers/counts/indptr'][:].astype(np.int64)
    payload=int(np.sum(ptr[selected+1]-ptr[selected]))*8
    return selected,np.isin(selected,rna),targets,payload


def read_selected_counts(node, columns, progress, max_bytes):
    """Read only selected CSC columns, with an explicit payload bound."""
    source_ptr=node['indptr'][:].astype(np.int64)
    lengths=source_ptr[columns+1]-source_ptr[columns]
    ptr=np.r_[0,np.cumsum(lengths)];total=int(ptr[-1])*8
    if total>max_bytes:
        raise ValueError('Predeclared sparse panel exceeds read budget')
    values=np.empty(int(ptr[-1]),dtype=np.float32)
    rows=np.empty(int(ptr[-1]),dtype=np.int32)
    digests={name:hashlib.sha256() for name in ['data','indices']}
    shape=tuple(map(int,node.attrs['shape']))
    for j,col in enumerate(columns):
        a,b=int(source_ptr[col]),int(source_ptr[col+1]);lo,hi=int(ptr[j]),int(ptr[j+1])
        values[lo:hi]=node['data'][a:b];rows[lo:hi]=node['indices'][a:b]
        v=values[lo:hi];r=rows[lo:hi]
        if not np.isfinite(v).all() or (v<=0).any() or not np.equal(v,np.floor(v)).all():
            raise ValueError('Selected count entries are not positive integers')
        if (r<0).any() or (r>=shape[0]).any() or np.any(np.diff(r)<=0):
            raise ValueError('Invalid selected row identity/order')
        digests['data'].update(memoryview(v).cast('B'));digests['indices'].update(memoryview(r).cast('B'))
        if j%20==0 or j+1==len(columns):progress('read_counts_panel',genes_done=j+1,genes_total=len(columns),logical_bytes=hi*8,logical_total=total)
    return values,rows,ptr,{k:v.hexdigest() for k,v in digests.items()}


def correlation(a, b):
    a, b = np.asarray(a,dtype=float), np.asarray(b,dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum()<3 or np.std(a[keep])==0 or np.std(b[keep])==0:
        return float('nan')
    return float(np.corrcoef(a[keep],b[keep])[0,1])


def response(vector, reference, method):
    if method == 'A_mean_log2FC':
        return np.log2((vector+1)/(reference+1))
    if method == 'B_mean_log1p_difference':
        return vector-reference
    raise ValueError('Unknown representation')


def compare_arms(entries, means, references, ref_sizes, method, n, seed):
    """Same input arms/features for shared- and independent-reference variants."""
    rows = entries[entries.arm.eq(0)].copy()
    records = []
    batches = sorted(references)
    for batch in batches:
        group = rows[rows.batch.eq(batch)]
        a, b = references[batch]
        size_a, size_b = ref_sizes[batch]
        shared = (a*size_a+b*size_b)/(size_a+size_b)
        eligible = group[group.role.eq('targeting')]
        control = group[group.role.eq('ntc')]
        rng = np.random.default_rng(stable_seed(seed,n,batch))
        def record(kind, left, right, target):
            x0,x1 = means[int(left.arm_id)],means[int(left.arm_id)+1]
            y0,y1 = means[int(right.arm_id)],means[int(right.arm_id)+1]
            base = dict(category=kind,batch=batch,arm_cells=n,method=method,target=target,
                        guide_a=left.guide,guide_b=right.guide,
                        expression_r=(correlation(x0,y1)+correlation(y0,x1))/2)
            for scheme,ra,rb in [('shared',shared,shared),('independent',a,b)]:
                first=correlation(response(x0,ra,method),response(y1,rb,method))
                if left.guide==right.guide:
                    value=first
                else:
                    second=correlation(response(y0,ra,method),response(x1,rb,method))
                    value=(first+second)/2
                records.append(dict(**base,reference_scheme=scheme,r=value))
        for _,row in group.iterrows():
            record('same_guide_targeting' if row.role=='targeting' else 'same_guide_ntc',row,row,row.target)
        for target, members in eligible.groupby('target',sort=True):
            for i,j in itertools.combinations(members.index,2):
                left,right=members.loc[i],members.loc[j]
                record('same_target_different_guides',left,right,target)
                other=eligible[~eligible.target.eq(target)]
                if len(other):
                    record('different_target_background',left,other.iloc[int(rng.integers(len(other)))],target)
                if len(control)>=2:
                    ix=rng.choice(len(control),2,replace=False)
                    record('ntc_pair_background',control.iloc[ix[0]],control.iloc[ix[1]],target)
    return pd.DataFrame(records)


def describe(values):
    values=np.asarray(values,dtype=float);values=values[np.isfinite(values)]
    if not len(values):return dict(n=0)
    return dict(n=len(values),median=float(np.median(values)),q25=float(np.quantile(values,.25)),
                q75=float(np.quantile(values,.75)),fraction_above_02=float(np.mean(values>.2)))


def run(config_path):
    config=json.loads(Path(config_path).read_text())
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();pid=os.getpid();source=Path(config['source'])
    lock=out/'owner.lock';lock.write_text(json.dumps(dict(pid=pid,uid=os.getuid())))
    def progress(stage,**kw):
        status=dict(stage=stage,pid=pid,uid=os.getuid(),elapsed_seconds=time.monotonic()-start,**kw)
        temp=out/'status.tmp';temp.write_text(json.dumps(status));temp.replace(out/'status.json')
        print(json.dumps(status),flush=True)
    def alarm(*_):raise TimeoutError('Declared public-data audit time budget reached')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(config.get('max_seconds',1800))
    try:
        before=source.stat()
        if (before.st_size,before.st_mtime_ns)!=(config['expected_bytes'],config['expected_mtime_ns']):
            raise ValueError('Source differs from accepted frozen identity')
        with h5py.File(source,'r') as h:
            progress('metadata')
            needed=['_index','gene_target','gRNA','batch','channel','n_gRNA','n_UMI_counts','n_genes',
                    'target_knockdown','exceeds_ntc_fdr','perturbation_cross_batch_significant']
            obs=pd.DataFrame({k:read_field(h['obs'][k]) for k in needed})
            for key in ['n_gRNA','n_UMI_counts','n_genes','target_knockdown']:
                obs[key]=pd.to_numeric(obs[key],errors='raise')
            if not obs['_index'].is_unique or not obs.n_gRNA.eq(1).all():
                raise ValueError('Ambiguous record identity or non-single-guide cohort')
            genes=np.asarray(read_field(h['var/_index']),dtype=str)
            gene_ids=np.asarray(read_field(h['var/gene_ids']),dtype=str)
            if not pd.Index(genes).is_unique or not pd.Index(gene_ids).is_unique:
                raise ValueError('Ambiguous RNA identity')
            full_table,full_map,ref_map,plans,ref_a,ref_b=make_plan(obs,config['arm_sizes'],config['seed'])
            full_table.to_csv(out/'guide_batch_coverage.csv.gz',index=False)
            reference_names=dict(a=ref_a,b=ref_b)
            (out/'reference_guides.json').write_text(json.dumps(reference_names,indent=2))
            source_shape=tuple(map(int,h['layers/counts'].attrs['shape']))
            if source_shape!=(len(obs),len(genes)):
                raise ValueError('RNA count axes differ from cell/gene identities')
            panel_mode=config.get('feature_mode')=='preexisting_hvg_panel'
            panel_feature=None;target_panel=None
            if panel_mode:
                columns,panel_feature,target_panel,payload=select_panel(h,obs,config['panel_size'],config['target_panel_size'])
                progress('planned_panel',genes=len(columns),correlation_panel=config['panel_size'],logical_payload_bytes=payload)
                data,indices,ptr,digests=read_selected_counts(h['layers/counts'],columns,progress,config['max_payload_bytes'])
                genes=genes[columns];gene_ids=gene_ids[columns]
                totals=obs.n_UMI_counts.to_numpy(dtype=float)
                detected=obs.n_genes.to_numpy(dtype=np.int64)
                if (totals<=0).any():raise ValueError('Invalid author full-RNA totals')
                selected_totals=np.bincount(indices,weights=data,minlength=len(obs))
                selected_detected=np.bincount(indices,minlength=len(obs))
                if np.any(selected_totals>totals) or np.any(selected_detected>detected):
                    raise ValueError('Selected integer counts exceed author whole-cell totals')
            else:
                data,indices,ptr,totals,detected,digests=read_counts(h['layers/counts'],progress)
        progress('aggregate_profiles')
        means_raw=np.array([float(data[ptr[j]:ptr[j+1]].sum(dtype=np.float64))/len(obs) for j in range(len(genes))])
        feature=means_raw>config.get('raw_mean_feature_threshold',.5)
        if panel_mode:feature &= panel_feature
        target_gene=np.isin(genes,target_panel if panel_mode else obs.loc[~obs.gene_target.eq('NTC'),'gene_target'].unique())
        selected=np.flatnonzero(feature|target_gene)
        mask=feature[selected]
        if mask.sum()<100:
            raise ValueError('Too few expressed RNA features for profile comparison')
        pd.DataFrame(dict(gene=genes,gene_id=gene_ids,raw_mean=means_raw,correlation_feature=feature,
                          target_gene=target_gene)).to_csv(out/'features.csv.gz',index=False)
        batches=sorted(obs.batch.unique());ntc=obs.gene_target.eq('NTC').to_numpy()
        scale={b:float(np.median(totals[ntc&obs.batch.eq(b).to_numpy()])) for b in batches}
        factors=np.array([scale[b] for b in obs.batch])/totals
        cpm_factor=10000/totals
        n_ref=2*len(batches);ref_n=np.bincount(ref_map[ref_map>=0],minlength=n_ref)
        full_n=np.bincount(full_map,minlength=len(full_table))
        full_linear=np.zeros((len(full_table),len(selected)),dtype=np.float32)
        ref_linear=np.zeros((n_ref,len(selected)),dtype=np.float32)
        ref_log=np.zeros_like(ref_linear)
        arrays={n:{name:np.zeros((len(p['entries']),len(selected)),dtype=np.float32)
                   for name in ['linear','log']} for n,p in plans.items()}
        mt=np.zeros(len(obs),dtype=np.float64)
        for j,gene in enumerate(genes):
            if gene.upper().startswith('MT-'):
                lo,hi=ptr[j:j+2];mt+=np.bincount(indices[lo:hi],weights=data[lo:hi],minlength=len(obs))
        for col,j in enumerate(selected):
            lo,hi=ptr[j:j+2];rr=indices[lo:hi];raw=data[lo:hi].astype(np.float64)
            linear=raw*factors[rr];logged=np.log1p(raw*cpm_factor[rr])
            full_linear[:,col]=aggregate_column(linear,rr,full_map,len(full_table))/full_n
            ref_linear[:,col]=aggregate_column(linear,rr,ref_map,n_ref)/ref_n
            ref_log[:,col]=aggregate_column(logged,rr,ref_map,n_ref)/ref_n
            for n,plan in plans.items():
                for name,v in [('linear',linear),('log',logged)]:
                    arrays[n][name][:,col]=aggregate_column(v,rr,plan['assignment'],len(plan['entries']))/n
            if col%250==0:progress('aggregate_profiles',genes_done=col,genes_total=len(selected))
        del data,indices
        progress('evaluate_correlations')
        tables=[]
        ref_sizes={b:(int(ref_n[bi*2]),int(ref_n[bi*2+1])) for bi,b in enumerate(batches)}
        for n,plan in plans.items():
            plan['entries'].to_csv(out/f'arms_n{n}.csv.gz',index=False)
            for name,key,references in [('A_mean_log2FC','linear',ref_linear),('B_mean_log1p_difference','log',ref_log)]:
                refs={b:(references[bi*2,mask],references[bi*2+1,mask]) for bi,b in enumerate(batches)}
                table=compare_arms(plan['entries'],arrays[n][key][:,mask],refs,ref_sizes,name,n,config['seed'])
                tables.append(table)
                np.savez_compressed(out/f'profiles_n{n}_{key}.npz',means=arrays[n][key],genes=genes[selected],
                                    correlation_feature=mask,reference_means=references)
        pairs=pd.concat(tables,ignore_index=True)
        pairs.to_csv(out/'correlation_pairs.csv.gz',index=False)
        stats=[]
        for keys,g in pairs.groupby(['method','arm_cells','reference_scheme','category'],sort=True):
            stats.append(dict(zip(['method','arm_cells','reference_scheme','category'],keys),**describe(g.r)))
        summary_table=pd.DataFrame(stats);summary_table.to_csv(out/'correlation_summary.csv',index=False)
        effect_rows=[];lookup={name:j for j,name in enumerate(genes[selected])}
        for row in full_table.itertuples():
            if row.target=='NTC' or row.target not in lookup or (panel_mode and row.target not in target_panel):continue
            col=lookup[row.target];bi=batches.index(row.batch)
            a,b=ref_n[bi*2:bi*2+2]
            reference=(ref_linear[bi*2,col]*a+ref_linear[bi*2+1,col]*b)/(a+b)
            value=float(full_linear[row.group,col])
            effect_rows.append(dict(batch=row.batch,guide=row.guide,target=row.target,cells=row.cells,
                                   target_normalized_mean=value,ntc_normalized_mean=float(reference),
                                   relative_change=value/reference-1 if reference>0 else np.nan,
                                   log2FC_pseudocount_1=float(np.log2((value+1)/(reference+1)))))
        effect_table=pd.DataFrame(effect_rows);effect_table.to_csv(out/'target_RNA_changes.csv.gz',index=False)
        qc=pd.DataFrame(dict(cell_barcode=obs['_index'],batch=obs.batch,guide=obs.gRNA,target=obs.gene_target,
                             total_counts=totals,detected_RNA=detected,mt_percent=mt/totals*100,
                             author_total_counts=obs.n_UMI_counts,author_detected_RNA=obs.n_genes))
        qc.to_csv(out/'cell_count_QC.csv.gz',index=False)
        proof=dict(source=str(source),source_bytes=int(before.st_size),source_mtime_ns=int(before.st_mtime_ns),
                   historical_container_sha256=config['historical_sha256'],fresh_whole_container_hash=False,
                   sparse_payload_sha256=digests,rows=len(obs),genes=len(genes),counts_layer='counts',
                   read_scope='selected_RNA_columns' if panel_mode else 'all_RNA_columns',
                   normalization_totals_source='author_obs' if panel_mode else 'recomputed_all_counts',
                   source_genes=source_shape[1],mitochondrial_counts_recomputed=True,
                   selected_gene_ids=gene_ids.tolist() if panel_mode else None,
                   target_RNA_panel=target_panel,
                   author_total_count_mismatches=None if panel_mode else int(np.count_nonzero(totals!=obs.n_UMI_counts.to_numpy())),
                   author_detected_RNA_mismatches=None if panel_mode else int(np.count_nonzero(detected!=obs.n_genes.to_numpy())),
                   source_unchanged=True)
        after=source.stat()
        if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
            raise ValueError('Source identity changed during audit')
        (out/'input_proof.json').write_text(json.dumps(proof,indent=2))
        effect_eligible=effect_table[effect_table.cells.ge(50)&effect_table.ntc_normalized_mean.ge(.5)]
        final=dict(scope='author_selected_strong_perturbation_cohort_only',feature_scope='fixed_preexisting_HVG_panel' if panel_mode else 'all_RNA_raw_mean_filter',config=config,
                   cells=len(obs),source_genes=source_shape[1],gene_columns_read=len(genes),ntc_cells=int(ntc.sum()),targets=int(obs.loc[~ntc,'gene_target'].nunique()),
                   correlation_features=int(mask.sum()),normalized_NTC_median_by_batch=scale,
                   qc_quantiles={k:qc[k].quantile([.25,.5,.75]).to_dict() for k in ['total_counts','detected_RNA','mt_percent']},
                   target_RNA_changes=describe(effect_eligible.relative_change),
                   fraction_target_RNA_lower=float(effect_eligible.relative_change.lt(0).mean()),
                   eligible_target_RNA_groups=len(effect_eligible),
                   author_flags={k:obs[k].astype(str).value_counts().to_dict() for k in ['exceeds_ntc_fdr','perturbation_cross_batch_significant']},
                   elapsed_seconds=time.monotonic()-start,
                   max_rss_kib=__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss,
                   code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip(),
                   versions={m:__import__(m).__version__ for m in ['numpy','pandas','h5py']},
                   inference='Descriptive retained-cohort reproducibility; no full-atlas QC pass rate, independent validation, DNA editing rate, or calibrated hit calls.')
        (out/'summary.json').write_text(json.dumps(final,indent=2))
        progress('complete')
    except BaseException as error:
        progress('failed',error=repr(error));raise
    finally:
        signal.alarm(0);lock.unlink(missing_ok=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True)
    run(parser.parse_args().config)
