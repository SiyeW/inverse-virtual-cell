"""Evaluate retention and RNA effects after the declared RNA-repression gates.

Reconstruct selected means by subtracting only removed cells from accepted full
guide means. Counts for each removed cell are checkpointed before evaluation.
No additional NTC cleaning, new DE significance gate, or DNA-efficiency claim.
"""
from pathlib import Path
from itertools import combinations
import argparse, hashlib, io, json, os, signal, time, zipfile
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from audit_knockdown_selection import screen
from compare_agreement_overlay import atomic, members, sha, apply_mask
from perturbseq_core_effects import pair_effects


def revised_mean(original, original_n, removed, remaining_n):
    original = np.asarray(original, float)
    removed = np.asarray(removed, float)
    if remaining_n <= 0 or removed.shape != original.shape or original_n < remaining_n:
        raise ValueError('Invalid subtraction axes or cell counts')
    value = (original * original_n - removed) / remaining_n
    if not np.isfinite(value).all() or value.min(initial=0) < -1e-8:
        raise ValueError('Removed contributions exceed accepted full means')
    return np.maximum(value, 0)


def reference_row(path, row):
    # Only the prefix through this row is consumed: do not claim full-member CRC.
    with zipfile.ZipFile(path) as z, z.open('reference_means.npy') as f:
        version = np.lib.format.read_magic(f)
        read_header = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
        shape, fortran, dtype = read_header(f)
        if fortran or len(shape) != 2 or dtype.kind != 'f':
            raise ValueError('Unexpected reference array layout')
        body = f.read((row+1)*shape[1]*dtype.itemsize)
        if len(body) != (row+1)*shape[1]*dtype.itemsize:
            raise ValueError('Truncated reference prefix')
        value = np.frombuffer(body, dtype=dtype).reshape(row+1, shape[1])[-1].copy()
    return value, dict(member='reference_means.npy', prefix_payload_bytes=len(body),
                       prefix_sha256=hashlib.sha256(body).hexdigest(), whole_member_crc_checked=False)


def source_identity(path):
    s = Path(path).stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev)


def strings(node):
    return node.asstr()[()]


def read_removed(config, cells, genes, requests, output, progress):
    """Use cell-major sparse arrays; never read the whole count payload."""
    source = Path(config['path'])
    if source_identity(source) != config['identity']:
        raise ValueError('Partial-read source identity changed')
    output.mkdir(exist_ok=True)
    results = {}
    if not requests:
        return results, dict(cells_read=0, raw_count_container_fully_read=False)
    requested_ids = genes.gene_id.astype(str).to_numpy()
    consumed = 0
    with h5py.File(source, 'r', rdcc_nbytes=16*1024**2) as f:
        if config['format'] == 'csr_h5ad_counts':
            node = f['layers']['counts']
            if node.attrs.get('encoding-type') != 'csr_matrix':
                raise ValueError('Expected cell-major CSR counts layer')
            barcodes = strings(f['obs'][f['obs'].attrs['_index']])
            all_genes = strings(f['var'][f['var'].attrs['_index']])
            rna_source = np.arange(len(all_genes))
        elif config['format'] == '10x_csc':
            node = f['matrix']
            barcodes = strings(node['barcodes'])
            all_genes = strings(node['features']['id'])
            rna_source = np.flatnonzero(strings(node['features']['feature_type']) == 'Gene Expression')
        else:
            raise ValueError('Unsupported sparse source')
        if not pd.Index(barcodes).is_unique or not pd.Index(all_genes).is_unique:
            raise ValueError('Ambiguous source axes')
        gene_positions = pd.Index(all_genes).get_indexer(requested_ids)
        if (gene_positions < 0).any() or not np.array_equal(np.sort(gene_positions), rna_source):
            raise ValueError('Frozen full RNA universe differs')
        remap = np.full(len(all_genes), -1, dtype=np.int32)
        remap[gene_positions] = np.arange(len(genes))
        source_rows = pd.Index(barcodes).get_indexer(cells.cell_barcode.astype(str))
        if (source_rows < 0).any():
            raise ValueError('Frozen cell barcode absent from count source')
        ptr = node['indptr'][()]
        data, indices = node['data'], node['indices']
        for i, q in enumerate(requests):
            row = q['row']; cell = cells.iloc[row]; bc = str(cell.cell_barcode)
            key = hashlib.sha256(bc.encode()).hexdigest()
            p = output/(key+'.npz'); receipt = p.with_suffix('.json')
            contract = dict(source=config, barcode=bc, row=row, total=int(cell.total_counts),
                            gene_axis_sha256=hashlib.sha256('\n'.join(requested_ids).encode()).hexdigest())
            if receipt.exists():
                r = json.loads(receipt.read_text())
                if r['contract'] != contract or sha(p) != r['sha256']:
                    raise ValueError('Removed-cell checkpoint changed')
                with np.load(p, allow_pickle=False) as z:
                    ix, counts = z['indices'], z['counts']
            else:
                a, b = map(int, ptr[source_rows[row]:source_rows[row]+2])
                counts = np.asarray(data[a:b]); source_ix = indices[a:b]
                consumed += counts.nbytes + source_ix.nbytes
                if consumed > 64*1024**2:
                    raise ValueError('Unexpectedly large removed-cell read')
                mapped = remap[source_ix]; use = mapped >= 0
                ix, counts = mapped[use], counts[use]
                if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
                    raise ValueError('Sparse source is not nonnegative integer counts')
                counts = counts.astype(np.int64)
                if int(counts.sum()) != int(cell.total_counts):
                    raise ValueError('Removed-cell full RNA total differs')
                # The compact accepted panel supplies a second, exact count check.
                hit = counts[ix == q['target_column']].sum()
                if int(hit) != q['target_count']:
                    raise ValueError('Removed-cell target UMI differs from accepted panel')
                with p.with_suffix('.tmp').open('wb') as f_out:
                    np.savez_compressed(f_out, indices=ix, counts=counts)
                p.with_suffix('.tmp').replace(p)
                atomic(receipt, dict(contract=contract, sha256=sha(p)))
            vector = sparse.csr_matrix((counts, ix, [0, len(ix)]), shape=(1, len(genes)))
            results[row] = vector
            if (i+1) % 32 == 0 or i+1 == len(requests):
                progress('removed_cell_counts', cells_done=i+1, cells_total=len(requests), logical_array_bytes=consumed)
    if source_identity(source) != config['identity']:
        raise ValueError('Count source changed during partial reads')
    return results, dict(cells_read=len(requests), logical_array_bytes=consumed,
                         source_identity=config['identity'], raw_count_container_fully_read=False,
                         per_cell_totals_and_target_UMIs_exact=True)


def effects(means, center, records, ntc_n):
    refs = np.broadcast_to(center, means.shape).copy()
    for i in np.flatnonzero(records.source_control_type.eq('non_targeting_control')):
        n = int(records.cells.iloc[i])
        refs[i] = (center*ntc_n - means[i]*n)/(ntc_n-n)
    return apply_mask(np.log2((means+1)/(refs+1)), records)


def median(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if len(values) else None


def pair_plan(records, selected):
    """Freeze same/other-target identities using initial cell counts, never RNA."""
    selected = records.loc[records.guide.isin(selected)].sort_values('guide')
    ordinary = selected.loc[selected.source_control_type.eq('targeting')]
    candidates = [(a, b) for a, b in combinations(ordinary.itertuples(), 2) if a.target != b.target]
    sizes = np.asarray([np.sort([a.cells,b.cells]) for a,b in candidates],float)
    ta = np.asarray([a.target for a,b in candidates],str)
    tb = np.asarray([b.target for a,b in candidates],str)
    rows = []
    for target, block in selected.groupby('target', sort=True):
        for a, b in combinations(block.itertuples(), 2):
            ordered = np.sort([a.cells, b.cells])
            allowed = (ta != target)&(tb != target)
            if allowed.any() and a.source_control_type == 'targeting':
                cost = np.abs(np.log(sizes/ordered)).sum(axis=1)
                cost[~allowed] = np.inf
                best = candidates[int(np.argmin(cost))]
                oa, ob = best[0].guide, best[1].guide
            else:
                oa, ob = '', ''
            rows.append(dict(target=target, guide_a=a.guide, guide_b=b.guide,
                             role=a.source_control_type, other_a=oa, other_b=ob))
    return pd.DataFrame(rows, columns=['target','guide_a','guide_b','role','other_a','other_b'])


def summarize_stage(label, profile, records, active, counts, center, plan, out):
    lookup = {g:i for i,g in enumerate(records.guide)}
    ntc = records.loc[records.source_control_type.eq('non_targeting_control') & records.cells.gt(25)]
    guide_rows = []
    for g in active:
        i = lookup[g]; v = profile[i]; good = np.isfinite(v)
        matched = ntc.assign(cost=np.abs(np.log(ntc.cells/int(records.cells.iloc[i])))).sort_values(['cost','guide']).head(4)
        comparison = [float(np.sqrt(np.mean(profile[lookup[n]][good]**2))) for n in matched.guide]
        rms = float(np.sqrt(np.mean(v[good]**2)))
        guide_rows.append(dict(stage=label, guide=g, target=records.target.iloc[i], cells=counts[g],
                               rms=rms, median_absolute_change=float(np.median(abs(v[good]))),
                               nearby_NTC_rms_median=median(comparison),
                               rms_above_nearby_NTC_median=rms > median(comparison)))
    gf = pd.DataFrame(guide_rows); gf.to_csv(out/(label+'_guide_effects.csv'), index=False)
    prs = []
    for row in plan.itertuples():
        a = pair_effects(profile[lookup[row.guide_a]], profile[lookup[row.guide_b]])
        other = pair_effects(profile[lookup[row.other_a]], profile[lookup[row.other_b]]) if row.other_a else {}
        prs.append(dict(stage=label, target=row.target, role=row.role, guide_a=row.guide_a, guide_b=row.guide_b,
                        same_r=a['pearson'], other_a=row.other_a, other_b=row.other_b,
                        other_r=other.get('pearson'), same_distance=a['mean_absolute_difference'],
                        other_distance=other.get('mean_absolute_difference')))
    pf = pd.DataFrame(prs); pf.to_csv(out/(label+'_pairs.csv'), index=False)
    ordinary = pf.loc[pf.role.eq('targeting')] if len(pf) else pf
    if len(ordinary):
        usable = ordinary.dropna(subset=['same_r','other_r'])
        agreement = dict(pairs=len(ordinary), targets=ordinary.target.nunique(),
                         same_r_median=median(ordinary.same_r), other_r_median=median(ordinary.other_r),
                         per_target_same_r_median=median(ordinary.groupby('target').same_r.median()),
                         paired_advantage_median=median(usable.same_r-usable.other_r),
                         same_r_higher_fraction=float((usable.same_r > usable.other_r).mean()) if len(usable) else None)
    else:
        agreement = dict(pairs=0, targets=0)
    return dict(guides=len(gf), rms_median=median(gf.rms), nearby_NTC_rms_median=median(gf.nearby_NTC_rms_median),
                rms_above_nearby_NTC_fraction=float(gf.rms_above_nearby_NTC_median.mean()), ordinary_agreement=agreement)


def dataset(d, out, progress):
    out.mkdir(exist_ok=True)
    root = Path(d['target_result'])
    for name, expected in d['cache_metadata_sha256'].items():
        if sha(root/name) != expected: raise ValueError('Saved count metadata changed')
    for name, key in [('target_counts.npz','target_counts_sha256'),('target_counts.json','target_counts_receipt_sha256'),('guide_target_RNA.csv','point_table_sha256')]:
        if sha(root/name) != d[key]: raise ValueError('Saved panel changed')
    cells = pd.read_csv(root/'cells.csv.gz', keep_default_na=False)
    genes = pd.read_csv(root/'genes.csv.gz', keep_default_na=False)
    if len(cells)!=d['expected_cells'] or int(cells.role.eq('non_targeting_control').sum())!=d['expected_NTC_cells'] or not cells.cell_barcode.is_unique:
        raise ValueError('Frozen main cohort differs')
    points = pd.read_csv(root/'guide_target_RNA.csv')
    receipt = json.loads((root/'target_counts.json').read_text())
    panel = sparse.load_npz(root/'target_counts.npz').tocsc()
    col = {str(genes.iloc[j].gene_id): k for k,j in enumerate(receipt['positions'])}
    gene_index = {g:i for i,g in enumerate(genes.gene_id.astype(str))}
    groups = cells.groupby('guide_id', sort=False).indices
    factors = 10000/cells.total_counts.to_numpy(float)
    selected = {}; requests = []; retained = []; stage = []
    for row in points.itertuples(index=False):
        use = groups.get(row.guide, np.array([], int)); target_col = col[row.target_gene_id]
        raw = panel[use, target_col].toarray().ravel()
        values = raw*factors[use]
        check = screen(values, row.NTC_mean_per_10000)
        if len(use) != int(row.cells) or (len(use) and not np.isclose(values.mean(), row.mean_per_10000, atol=1e-11, rtol=1e-11)):
            raise ValueError('Original target profile bridge failed')
        stage.append(dict(guide=row.guide, target=row.target, target_gene_id=row.target_gene_id,
                          cells=int(row.cells), NTC_detected_fraction=row.NTC_detected_fraction,
                          detected_cells=int(row.detected_cells), **check))
        if check['final_pass']:
            good = values <= row.NTC_mean_per_10000
            selected[row.guide] = int(good.sum())
            for k in use[good]: retained.append(dict(guide=row.guide, cell_barcode=str(cells.cell_barcode.iloc[k]), analysis_row=int(k)))
            for k, count in zip(use[~good], raw[~good]):
                requests.append(dict(row=int(k), guide=row.guide, target_column=gene_index[row.target_gene_id], target_count=int(count)))
    sf = pd.DataFrame(stage); sf.to_csv(out/'retention_by_guide.csv', index=False)
    pd.DataFrame(retained).to_csv(out/'retained_cells.csv.gz', index=False)
    cache = out/'baseline_profiles.npz'; cr = cache.with_suffix('.json')
    baseline = Path(d['baseline'])
    for p,h in d['metadata_sha256'].items():
        if sha(Path(p)) != h: raise ValueError('Original metadata changed')
    records = pd.read_csv(baseline/'tables/guide_effect_metadata.csv')
    if cr.exists():
        if sha(cache) != json.loads(cr.read_text())['sha256']: raise ValueError('Profile cache changed')
        with np.load(cache, allow_pickle=False) as z: old = {k:z[k] for k in z.files}
    else:
        old, proof = members(baseline/'tables/rna_effect_profiles.npz', ['means','guides','gene_ids'])
        if not np.array_equal(old['guides'], records.guide): raise ValueError('Profile row identities differ')
        ri = int(np.flatnonzero(records.source_control_type.ne('non_targeting_control'))[0])
        old['center'], refproof = reference_row(baseline/'tables/rna_effect_profiles.npz', ri)
        np.savez_compressed(cache, **old); atomic(cr, dict(sha256=sha(cache), source=proof, reference_prefix=refproof))
        progress('saved_original_means', compressed_bytes=sum(v['compressed_bytes'] for v in proof['consumed_members']))
    if not np.array_equal(old['guides'], records.guide): raise ValueError('Cached profile rows differ')
    old_effect = effects(old['means'], old['center'], records, d['expected_NTC_cells'])
    original_pairs = pd.read_csv(baseline/'tables/within_target_pairs.csv')
    lookup = {g:i for i,g in enumerate(records.guide)}
    bridge = 0.0
    for q in original_pairs.itertuples():
        r = pair_effects(old_effect[lookup[q.guide_a]],old_effect[lookup[q.guide_b]])['pearson']
        if pd.notna(q.pearson) and r is not None: bridge=max(bridge,abs(r-q.pearson))
    if bridge > 1e-10: raise ValueError('Original pair correlations did not reproduce')
    raw_rows, extraction = read_removed(d['sparse_source'], cells, genes, requests, out/'removed_cell_cache', progress)
    selected_columns = pd.Index(genes.gene_id.astype(str)).get_indexer(old['gene_ids'].astype(str))
    if (selected_columns<0).any(): raise ValueError('RNA profile features absent from count source')
    removals = {g:np.zeros(len(selected_columns)) for g in selected}
    for q in requests:
        row=q['row'];removals[q['guide']]+=raw_rows[row][:,selected_columns].toarray().ravel()*(d['baseline_scale']/cells.total_counts.iloc[row])
    new_means=old['means'].copy()
    for g,n in selected.items():
        i=lookup[g];new_means[i]=revised_mean(old['means'][i],int(records.cells.iloc[i]),removals[g],n)
    new_effect=effects(new_means,old['center'],records,d['expected_NTC_cells'])
    np.savez_compressed(out/'selected_profiles.npz',guides=old['guides'],gene_ids=old['gene_ids'],
                        original_means=old['means'],selected_means=new_means,center=old['center'])
    plan=pair_plan(records,set(selected));plan.to_csv(out/'frozen_pairs.csv',index=False)
    original_counts=dict(zip(records.guide,records.cells))
    stages={}
    for label,profile,active,counts in [('selected_guides_original_cells',old_effect,list(selected),original_counts),('selected_cells',new_effect,list(selected),selected)]:
        stages[label]=summarize_stage(label,profile,records,active,counts,old['center'],plan,out)
    broad=set(sf.loc[sf.cells.gt(25)&sf.reference_defined,'guide']);broad_plan=pair_plan(records,broad)
    stages['all_covered']=summarize_stage('all_covered',old_effect,records,sorted(broad),original_counts,old['center'],broad_plan,out)
    measured=sf.final_pass&sf.NTC_detected_fraction.ge(.1)
    measured_set=set(sf.loc[measured,'guide']);measured_plan=plan.loc[plan.guide_a.isin(measured_set)&plan.guide_b.isin(measured_set)]
    stages['selected_cells_detection_ge10']=summarize_stage('selected_cells_detection_ge10',new_effect,records,sorted(measured_set),selected,old['center'],measured_plan,out)
    all_cells=int(points.cells.sum());all_obs=int(points.cells.gt(0).sum());all_targets=points.loc[points.cells.gt(0),'target_gene_id'].nunique()
    retention=dict(designed_guides=len(points),designed_targets=points.target_gene_id.nunique(),
                   observed_guides=all_obs,observed_targets=int(all_targets),targeting_cells=all_cells,
                   admitted_guides=len(selected),admitted_targets=sf.loc[sf.final_pass,'target_gene_id'].nunique(),
                   retained_targeting_cells=sum(selected.values()),
                   guide_fraction_of_design=len(selected)/len(points),guide_fraction_of_observed=len(selected)/all_obs,
                   target_fraction_of_design=sf.loc[sf.final_pass,'target_gene_id'].nunique()/points.target_gene_id.nunique(),
                   cell_fraction_of_targeting=sum(selected.values())/all_cells,
                   removed_from_retained_guides=len(requests), NTC_cells_unchanged=d['expected_NTC_cells'])
    strata=[]
    for lo,hi in [(0,.01),(.01,.05),(.05,.1),(.1,1.01)]:
        sub=sf.loc[sf.final_pass&sf.NTC_detected_fraction.ge(lo)&sf.NTC_detected_fraction.lt(hi)]
        strata.append(dict(NTC_detection_min=lo,NTC_detection_max=hi,guides=len(sub),targets=sub.target_gene_id.nunique(),cells=int(sub.after_cell_gate.sum())))
    result=dict(dataset=d['id'],retention=retention,stages=stages,NTC_detection_strata=strata,
                original_pair_bridge_max_abs_error=bridge,extraction=extraction,
                original_RNA_features=len(old['gene_ids']),original_normalization_depth=d['baseline_scale'],
                median_RMS_effect_vector_change=median([np.sqrt(np.nanmean((new_effect[lookup[g]]-old_effect[lookup[g]])**2)) for g in selected]),
                RNA_selection_unchanged=True,reference_unchanged=True,paper_full_pipeline=False,
                no_DE_significance_calls=True, no_DNA_efficiency_claim=True)
    atomic(out/'summary.json',result)
    return result


def run(config):
    out=Path(config['output']);out.mkdir(exist_ok=True)
    if (out/'owner.lock').exists(): raise ValueError('Another owner or unreconciled prior writer')
    lock=out/'owner.lock';lock.write_text(json.dumps(dict(pid=os.getpid(),uid=os.getuid())))
    started=time.time()
    def progress(stage,**kw): atomic(out/'status.json',dict(stage=stage,elapsed_seconds=time.time()-started,**kw))
    def timeout(*_): raise TimeoutError('Bounded effect audit stopped; saved profiles and cells remain reusable')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(config.get('max_seconds',900))
    try:
        results=[]
        for d in config['datasets']:
            def update(stage,**kw): progress(stage,dataset=d['id'],**kw)
            results.append(dataset(d,out/d['id'],update))
        atomic(out/'summary.json',dict(commit=config['expected_commit'],datasets=results,elapsed_seconds=time.time()-started))
        progress('complete')
    except Exception as ex:
        progress('failed',error=repr(ex));raise
    finally:
        signal.alarm(0);lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');run(json.loads(Path(p.parse_args().config).read_text()))
