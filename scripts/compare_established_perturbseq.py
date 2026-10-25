"""Finite descriptive comparison of established Perturb-seq representations.

Uses Scanpy normalization/HVG/PCA/DE and sklearn nearest controls. This is a
method sensitivity study, not biological replicate inference or Mixscape KO
classification. Inputs and outputs are specified in an external run contract.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import io
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
from scipy import sparse


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def canonical_counts(x):
    x = x.tocsr(copy=False)
    if not x.has_sorted_indices or not x.has_canonical_format:
        raise ValueError('Canonical CSR required, not silently repaired')
    h = hashlib.sha256(np.asarray(x.shape, dtype='<i8').tobytes())
    for name, dtype in [('indptr', '<i8'), ('indices', '<i4'), ('data', '<i4')]:
        values = getattr(x, name)
        h.update(name.encode())
        for start in range(0, len(values), 1048576):
            v = values[start:start + 1048576]
            if name == 'data' and (not np.isfinite(v).all() or (v < 0).any() or
                                   not np.equal(v, np.floor(v)).all()):
                raise ValueError('Noninteger or negative counts')
            h.update(v.astype(dtype, copy=False).tobytes())
    return h.hexdigest()


def identity(path):
    s = Path(path).stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns,
                device=s.st_dev, inode=s.st_ino)


def text_array(dataset):
    return dataset.asstr()[()]


def checked_index(available, requested):
    index = pd.Index(available)
    if not index.is_unique or not pd.Index(requested).is_unique:
        raise ValueError('Duplicate input identities')
    result = index.get_indexer(requested)
    if (result < 0).any():
        raise ValueError('Missing requested identity')
    return result


def load_counts(config, progress):
    """Use either count-identical native cache or hash-verified 10x container."""
    source = Path(config['source'])
    before = identity(source)
    cells = pd.read_csv(config['cells'], dtype={'cell_barcode': str, 'guide_id': str})
    catalogue = pd.read_csv(config['catalogue'], keep_default_na=False)
    if not cells.cell_barcode.is_unique or not catalogue.guide.is_unique:
        raise ValueError('Ambiguous cell/catalogue identity')
    if config['input_format'] == 'transposed_scanpy_cache':
        contract = json.loads(Path(config['contract']).read_text())
        for name, entry in contract['files'].items():
            if sha_file(Path(config['contract']).parent / name) != entry['sha256']:
                raise ValueError('Input metadata hash differs')
        genes = pd.read_csv(config['genes'], keep_default_na=False)
        features = pd.read_csv(config['source_features'], sep='\t', header=None,
                               names=['gene_id', 'gene_name', 'feature_type'])
        with h5py.File(source, 'r') as f:
            # scanpy.read_mtx caches the original genes-by-cells matrix.
            node = f['X']
            shape = tuple(node.attrs['shape'])
            raw = sparse.csr_matrix((node['data'][()], node['indices'][()],
                                     node['indptr'][()]), shape=shape)
            cache_barcodes = text_array(f['var']['_index'])
            cache_names = text_array(f['obs']['_index'])
        if len(features) != shape[0]:
            raise ValueError('Cache/features cardinality differs')
        # Matrix Market reader uses positional row names before gene annotation.
        if not np.array_equal(cache_names, np.arange(shape[0]).astype(str)):
            raise ValueError('Unexpected cache row identity')
        rows = checked_index(cache_barcodes, cells.cell_barcode)
        cols = checked_index(features.gene_id, genes.gene_id)
        x = raw[cols][:, rows].T.tocsr()
        del raw
        x.sort_indices()
        got = canonical_counts(x)
        if got != contract['canonical_counts_sha256']:
            raise ValueError('Native cache counts differ from authoritative cohort')
        proof = dict(canonical_counts_sha256=got, authoritative_contract=contract,
                     cache_count_identity_verified=True)
    elif config['input_format'] == '10x_h5':
        # One bounded sequential container read avoids repeated small NAS reads.
        payload = bytearray()
        h = hashlib.sha256()
        with source.open('rb') as stream:
            while True:
                block = stream.read(4 * 1024**2)
                if not block:
                    break
                payload.extend(block)
                h.update(block)
                if len(payload) // (64 * 1024**2) != (len(payload)-len(block)) // (64 * 1024**2):
                    progress('read_input', bytes=len(payload), total=before['bytes'])
        if h.hexdigest() != config['source_sha256']:
            raise ValueError('10x source container SHA256 mismatch')
        with h5py.File(io.BytesIO(payload), 'r') as f:
            node = f['matrix']
            genes = pd.DataFrame({k: text_array(node['features'][v]) for k, v in
                                  [('gene_id', 'id'), ('gene_name', 'name'),
                                   ('feature_type', 'feature_type')]})
            all_barcodes = text_array(node['barcodes'])
            rows = checked_index(all_barcodes, cells.cell_barcode)
            shape = tuple(node['shape'][()])
            data = node['data'][()]
            if (data < 0).any() or data.max() > np.iinfo(np.int32).max:
                raise ValueError('Counts outside supported range')
            data = data.astype(np.int32)
            raw = sparse.csc_matrix((data, node['indices'][()], node['indptr'][()]), shape=shape)
            x = raw[:, rows].T.tocsr()
        del raw, payload, data
        rna = genes.feature_type.eq('Gene Expression').to_numpy()
        x = x[:, rna].tocsr()
        genes = genes.loc[rna].reset_index(drop=True)
        proof = dict(source_sha256=h.hexdigest(), whole_source_hash_verified=True)
    else:
        raise ValueError('Unsupported declared input')
    if identity(source) != before:
        raise ValueError('Source identity changed')
    totals = np.asarray(x.sum(1)).ravel()
    if not np.array_equal(totals, cells.total_counts.to_numpy()):
        raise ValueError('Authoritative per-cell RNA totals differ')
    if len(cells) != config['expected_cells'] or x.shape[1] != config['expected_genes']:
        raise ValueError('Cohort shape differs')
    annotation = catalogue.set_index('guide').reindex(cells.guide_id)
    if annotation.source_control_type.isna().any():
        raise ValueError('Missing guide design annotation')
    cells['target'] = annotation.target_gene_name.to_numpy()
    cells['role'] = annotation.source_control_type.to_numpy()
    # Annotations supply design identity only, never legacy coverage/selection.
    proof.update(source=str(source), source_identity=before,
                 source_modified=False, cells=len(cells), genes=x.shape[1], nnz=x.nnz)
    return x, cells, genes, proof


def pearson(a, b, mask=None):
    keep = np.isfinite(a) & np.isfinite(b)
    if mask is not None:
        keep &= mask
    a, b = np.asarray(a)[keep], np.asarray(b)[keep]
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def ntc_z_parameters(normalized, ntc):
    controls = normalized[ntc]
    mean = np.asarray(controls.mean(0)).ravel().astype(float)
    second = np.asarray(controls.multiply(controls).mean(0)).ravel().astype(float)
    n = controls.shape[0]
    variance = np.maximum(0, second - mean**2) * n / (n-1)
    sd = np.sqrt(variance)
    return mean, sd, sd > 1e-6


def nearest_controls(embedding, controls, guide_labels, k=20):
    """Exact documented local-NTC primitive, excluding whole own NTC guide."""
    from sklearn.neighbors import NearestNeighbors
    controls = np.asarray(controls, bool)
    ix = np.flatnonzero(controls)
    guide_labels = np.asarray(guide_labels)
    max_own = int(pd.Series(guide_labels[ix]).value_counts().max())
    if len(ix)-max_own < k:
        raise ValueError('Insufficient disjoint NTC neighbors')
    model = NearestNeighbors(n_neighbors=k+max_own, algorithm='brute', n_jobs=2).fit(embedding[ix])
    result = np.empty((len(embedding), k), dtype=int)
    for start in range(0, len(embedding), 512):
        query = model.kneighbors(embedding[start:start+512], return_distance=False)
        for j, candidate in enumerate(query):
            row = start+j
            candidate = ix[candidate]
            eligible = candidate[guide_labels[candidate] != guide_labels[row]]
            result[row] = eligible[:k]
    return result


def profiles(matrix, groups, control_mask, guides, *, mode='difference', sd=None):
    """Full group leave-out references; never include an NTC group's own cells."""
    records, output = [], []
    control_sum = np.asarray(matrix[control_mask].sum(0)).ravel().astype(float)
    n_controls = int(control_mask.sum())
    for name, rows in groups:
        rows = np.asarray(rows)
        own_guides = np.unique(guides[rows])
        removed = control_mask & np.isin(guides, own_guides)
        n_ref = n_controls - int(removed.sum())
        if n_ref < 2:
            raise ValueError('Reference empty after group exclusion')
        a = np.asarray(matrix[rows].mean(0)).ravel()
        removed_sum = np.asarray(matrix[removed].sum(0)).ravel() if removed.any() else 0
        b = (control_sum - removed_sum) / n_ref
        value = np.log2((a+1)/(b+1)) if mode == 'lfc' else a-b
        if sd is not None:
            value = np.divide(value, sd, out=np.full_like(value, np.nan), where=sd>1e-6)
        records.append(dict(group=name, cells=len(rows), reference_cells=n_ref))
        output.append(value)
    return pd.DataFrame(records), np.asarray(output)


def summarize(values):
    a = np.asarray(values, float); a = a[np.isfinite(a)]
    return dict(n=len(a), median=float(np.median(a)) if len(a) else None,
                p10=float(np.quantile(a, .1)) if len(a) else None,
                p90=float(np.quantile(a, .9)) if len(a) else None)


def group_sets(cells, seed):
    rng = np.random.default_rng(seed)
    guide_groups = [(name, np.asarray(rows)) for name, rows in
                    cells.groupby('guide_id', sort=True).indices.items() if len(rows)>=20]
    half_groups = []
    for name, rows in guide_groups:
        order = rng.permutation(rows)
        half_groups.extend([(name+'::A', order[::2]), (name+'::B', order[1::2])])
    # Target pooling is separate from guide replication; controls stay separate.
    target_groups = [(name, np.asarray(rows)) for name, rows in
                     cells.loc[cells.role.ne('non_targeting_control')].groupby('target', sort=True).groups.items()
                     if name and len(rows)>=40]
    return guide_groups, half_groups, target_groups


def evaluate(method, matrix, cells, genes, groups, mode, sd, feature_ids):
    ntc = cells.role.eq('non_targeting_control').to_numpy()
    guides = cells.guide_id.to_numpy()
    full_groups, half_groups, target_groups = groups
    frames = []
    for level, grouped in [('guide', full_groups), ('half', half_groups), ('target', target_groups)]:
        table, effects = profiles(matrix, grouped, ntc, guides, mode=mode, sd=sd)
        if mode == 'local':
            # matrix already consists of individual local-control residuals.
            effects = np.asarray([np.asarray(matrix[rows].mean(0)).ravel() for _, rows in grouped])
        table['level'] = level
        frames.append((table, effects))
    full, effect = frames[0]
    meta = cells.drop_duplicates('guide_id').set_index('guide_id')
    full['target'] = meta.loc[full.group, 'target'].to_numpy()
    full['role'] = meta.loc[full.group, 'role'].to_numpy()
    same, unrelated, ntc_pairs = [], [], []
    masks = {}
    for target in full.target.unique():
        # Consistent with old plot: target's own RNA is excluded from trans-response.
        masks[target] = genes.gene_name.to_numpy() != target
    candidates = []
    for i,j in itertools.combinations(range(len(full)), 2):
        a,b = full.iloc[i], full.iloc[j]
        mask = masks[a.target] & masks[b.target]
        if a.role == b.role == 'non_targeting_control':
            ntc_pairs.append(dict(a=a.group,b=b.group,r=pearson(effect[i],effect[j],mask)))
        elif a.role != 'non_targeting_control' and b.role != 'non_targeting_control':
            if a.target == b.target:
                same.append(dict(a=a.group,b=b.group,target=a.target,role=a.role,r=pearson(effect[i],effect[j],mask),n_a=int(a.cells),n_b=int(b.cells)))
            else:
                candidates.append((i,j))
    # Fixed seeded candidate pool; match log cell counts without looking at effects.
    rng = np.random.default_rng(731)
    candidates = np.asarray(candidates)
    if len(candidates):
        selected = rng.choice(len(candidates), min(4000,len(candidates)), replace=False)
        pool=candidates[selected]
        pool_sizes=np.sort(np.log(full.cells.to_numpy()[pool]),axis=1)
        selected_matches=[]
        for pair in same:
            wanted=np.sort(np.log([pair['n_a'],pair['n_b']]))
            best=int(np.argmin(np.square(pool_sizes-wanted).sum(1)))
            selected_matches.append((pool[best],pair))
        for (i,j),matched in selected_matches:
            a,b = full.iloc[i],full.iloc[j]
            unrelated.append(dict(a=a.group,b=b.group,r=pearson(effect[i],effect[j],masks[a.target]&masks[b.target]),n_a=int(a.cells),n_b=int(b.cells),matched_to_a=matched['a'],matched_to_b=matched['b']))
    half, half_effect = frames[1]
    repeat = []
    for i in range(0,len(half),2):
        name = half.iloc[i].group.removesuffix('::A')
        row = meta.loc[name]
        repeat.append(dict(guide=name,target=row.target,role=row.role,cells=int(full.loc[full.group.eq(name),'cells'].iloc[0]),
                           r=pearson(half_effect[i],half_effect[i+1],masks[row.target])))
    target_table,target_effect=frames[2]
    target_pairs=[]
    # Disjoint guide sets per target, not cell splits of a mixture of identical guides.
    for target in target_table.group:
        eligible = full.loc[full.target.eq(target)&full.role.ne('non_targeting_control')]
        if len(eligible)<2:
            continue
        names = sorted(eligible.group)
        first = cells.guide_id.isin(names[::2]).to_numpy()
        second = cells.guide_id.isin(names[1::2]).to_numpy()
        split = [(target+'::guidesA',np.flatnonzero(first)),(target+'::guidesB',np.flatnonzero(second))]
        _, pair = profiles(matrix,split,ntc,guides,mode=mode,sd=sd)
        if mode=='local':
            pair = np.asarray([np.asarray(matrix[rows].mean(0)).ravel() for _,rows in split])
        target_pairs.append(dict(target=target,n_a=int(first.sum()),n_b=int(second.sum()),
                                 guides_a=len(names[::2]),guides_b=len(names[1::2]),
                                 r=pearson(pair[0],pair[1],masks[target])))
    tables=dict(same_target_pairs=pd.DataFrame(same),unrelated_pairs=pd.DataFrame(unrelated),
                ntc_pairs=pd.DataFrame(ntc_pairs),same_guide_halves=pd.DataFrame(repeat),
                target_disjoint_guides=pd.DataFrame(target_pairs),guide_groups=full)
    summary={key:summarize(table.r) for key,table in tables.items() if 'r' in table}
    summary.update(method=method,features=len(feature_ids),eligible_guides=len(full))
    ordinary=tables['same_target_pairs']
    summary['ordinary_same_target']=summarize(ordinary.loc[ordinary.role.eq('targeting'),'r']) if same else summarize([])
    summary['same_target_median_of_target_medians']=float(tables['same_target_pairs'].groupby('target').r.median().median()) if same else None
    return summary,tables


def run(config_path):
    import scanpy as sc
    import anndata as ad
    from sklearn.decomposition import PCA
    from threadpoolctl import threadpool_limits

    config=json.loads(Path(config_path).read_text())
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    lock=out/'owner.lock';lock.write_text(str(os.getpid()))
    start=time.monotonic()
    def progress(stage,**kwargs):
        status=dict(stage=stage,pid=os.getpid(),elapsed_seconds=time.monotonic()-start,**kwargs)
        tmp=out/'status.tmp';tmp.write_text(json.dumps(status));tmp.replace(out/'status.json')
        print(json.dumps(status),flush=True)
    def timeout(signum,frame):
        raise TimeoutError('Declared analysis wall-time budget reached')
    signal.signal(signal.SIGALRM,timeout)
    signal.alarm(config.get('max_seconds',900))
    try:
        with threadpool_limits(limits=2):
            progress('load_counts')
            x,cells,genes,proof=load_counts(config,progress)
            (out/'input_proof.json').write_text(json.dumps(proof,indent=2))
            ntc=cells.role.eq('non_targeting_control').to_numpy()
            if not ntc.any():raise ValueError('No NTCs')
            totals=np.asarray(x.sum(1)).ravel()
            scale=float(np.median(totals[ntc]))
            old_features=np.flatnonzero(np.asarray(x.mean(0)).ravel()>.5)
            norm=x.astype(np.float32,copy=False)
            del x
            norm.data*=np.repeat((scale/totals).astype(np.float32),np.diff(norm.indptr))
            groups=group_sets(cells,20261026)
            old=norm[:,old_features].tocsr()
            # Z-standardization uses the same linear normalization as the paper.
            progress('control_variance')
            control_mean,control_sd,valid_sd=ntc_z_parameters(norm,ntc)
            # log-normalize per cell before averaging, as in Scanpy peer workflow.
            norm.data=np.log1p(norm.data*(10000/scale)).astype(np.float32)
            a=ad.AnnData(norm,var=pd.DataFrame(index=genes.gene_id))
            sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=2000,inplace=True)
            hvg=np.flatnonzero(a.var.highly_variable.to_numpy()&valid_sd)
            hvg_log=norm[:,hvg].toarray()
            log_old=norm[:,old_features].tocsr()
            # Invert only selected features, avoiding a second full sparse copy.
            hvg_linear=np.expm1(hvg_log)*(scale/10000)
            del a,norm
            gc.collect()
            pd.DataFrame(dict(gene_id=genes.gene_id,old_mean_filter=np.isin(np.arange(len(genes)),old_features),
                              hvg=np.isin(np.arange(len(genes)),hvg),ntc_sd=control_sd)).to_csv(out/'features.csv.gz',index=False)
            progress('local_ntc_signature',hvg=len(hvg),old_features=len(old_features))
            embedding=PCA(n_components=40,svd_solver='randomized',random_state=20261026).fit_transform(hvg_log)
            neighbors=nearest_controls(embedding,ntc,cells.guide_id.to_numpy(),k=20)
            residual=hvg_log.copy()
            for first in range(0,len(cells),256):
                residual[first:first+256]-=hvg_log[neighbors[first:first+256]].mean(1)
            methods=[('A_original_mean_lfc',old,old_features,'lfc',None),
                     ('B_logmean_same_features',log_old,old_features,'difference',None),
                     ('C_logmean_hvg2000',hvg_log,hvg,'difference',None),
                     ('D_ntc_z_hvg2000',hvg_linear,hvg,'difference',control_sd[hvg]),
                     ('E_local_ntc_log_hvg2000',residual,hvg,'local',None)]
            summaries=[]
            for name,matrix,features,mode,sd in methods:
                progress('evaluate',method=name)
                summary,tables=evaluate(name,matrix,cells,genes.iloc[features],groups,mode,sd,features)
                summaries.append(summary)
                for key,table in tables.items():table.to_csv(out/(name+'_'+key+'.csv.gz'),index=False)
            # Conventional target-wise Wilcoxon pilot on fixed, label-blind HVGs.
            # NT controls are held-out whole guides; cell p-values are nominal only.
            progress('scanpy_target_de')
            coverage=cells.loc[~ntc].groupby('target').size().sort_values(ascending=False)
            positive=cells.loc[cells.role.eq('positive_control')].groupby('target').size().sort_values(ascending=False)
            chosen=list(dict.fromkeys(list(positive.index[:6])+list(coverage.index[:6])))
            negative=list(cells.loc[ntc].guide_id.value_counts().index[:6])
            labels=np.full(len(cells),'discard',dtype=object)
            labels[ntc]='Control'
            for name in chosen:labels[(~ntc)&cells.target.eq(name).to_numpy()]=name
            for name in negative:labels[cells.guide_id.eq(name).to_numpy()]='NTC::'+name
            take=labels!='discard'
            de=ad.AnnData(hvg_log[take],obs=pd.DataFrame({'group':pd.Categorical(labels[take])}),var=pd.DataFrame(index=genes.gene_id.iloc[hvg]))
            sc.tl.rank_genes_groups(de,groupby='group',reference='Control',method='wilcoxon',tie_correct=True,use_raw=False)
            result=sc.get.rank_genes_groups_df(de,group=None)
            result.to_csv(out/'scanpy_wilcoxon_pilot.csv.gz',index=False)
            stats=[]
            for name,table in result.groupby('group',observed=True):
                stats.append(dict(group=name,cells=int((labels==name).sum()),reference_cells=int((labels=='Control').sum()),
                                  nominal_bh05=int(table.pvals_adj.lt(.05).sum()),
                                  nominal_bh05_abs_lfc025=int((table.pvals_adj.lt(.05)&table.logfoldchanges.abs().gt(.25)).sum())))
            pd.DataFrame(stats).to_csv(out/'scanpy_wilcoxon_pilot_summary.csv',index=False)
            final=dict(config=config,representations=summaries,normalization_target_umi=scale,
                       genes=len(genes),cells=len(cells),ntc_cells=int(ntc.sum()),
                       input_scope=config['input_scope'],max_rss_kib=__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss,
                       elapsed_seconds=time.monotonic()-start,code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip(),
                       versions={p:__import__(p).__version__ for p in ['scanpy','anndata','numpy','scipy','sklearn']},
                       inference='descriptive and nominal per-cell pilot only; no biological replicate FDR or validated hits')
            (out/'summary.json').write_text(json.dumps(final,indent=2))
            progress('complete')
    except BaseException as error:
        progress('failed',error=repr(error));raise
    finally:
        signal.alarm(0);lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True)
    run(parser.parse_args().config)
