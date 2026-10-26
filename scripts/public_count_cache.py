"""Durable, axis-bound streaming cache of raw CSC counts (never dense X)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from perturbseq_core_io import atomic_json, sha256


def axis_hash(values):
    return hashlib.sha256(json.dumps(list(map(str, values)), ensure_ascii=True,
                                    separators=(',', ':')).encode()).hexdigest()


def column_spans(ptr, max_entries):
    """Contiguous complete-gene spans; a single large gene is never split."""
    ptr = np.asarray(ptr, dtype=np.int64)
    if len(ptr) < 2 or ptr[0] != 0 or (np.diff(ptr) < 0).any() or max_entries < 1:
        raise ValueError('Invalid CSC pointer/budget')
    first = 0
    while first < len(ptr)-1:
        end = int(np.searchsorted(ptr, int(ptr[first])+max_entries, side='right')-1)
        end = min(len(ptr)-1, max(first+1, end))
        yield first, end
        first = end


def check_arrays(data, rows, ptr, nrows):
    if len(data) != len(rows) or ptr[0] != 0 or ptr[-1] != len(data) or (np.diff(ptr) < 0).any():
        raise ValueError('Invalid sparse block shape')
    if not np.isfinite(data).all() or (data <= 0).any() or (data != np.floor(data)).any():
        raise ValueError('Raw positive integer counts required')
    if (data > np.iinfo(np.int32).max).any() or (rows < 0).any() or (rows >= nrows).any():
        raise ValueError('Count/row outside exact int32 bounds')
    # Canonical sorted unique row indices within every complete source column.
    if len(rows) > 1:
        boundary = np.zeros(len(rows)-1, dtype=bool)
        inside = ptr[1:-1]
        inside = inside[(inside > 0) & (inside < len(rows))]
        boundary[inside-1] = True
        if ((np.diff(rows) <= 0) & ~boundary).any():
            raise ValueError('Noncanonical or duplicate source cell entries')


def block_arrays(root, entry):
    path = Path(root)/entry['file']
    if path.stat().st_size != entry['bytes'] or sha256(path) != entry['sha256']:
        raise ValueError('Count checkpoint content changed: '+entry['file'])
    with np.load(path, allow_pickle=False) as block:
        return block['data'], block['rows'], block['ptr']


def cache_counts(node, root, source_key, row_ids, gene_ids, progress,
                 max_entries=8*1024**2):
    """An accepted block is verified/reused without rereading source payload.

    The caller exclusively owns root. On interrupted publication, unaccepted
    files are preserved and fail closed; they are never mistaken for receipts.
    """
    root = Path(root)
    root.mkdir(exist_ok=True)
    shape = tuple(map(int, node.attrs['shape']))
    encoding = node.attrs['encoding-type']
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding != 'csc_matrix' or shape != (len(row_ids), len(gene_ids)):
        raise ValueError('Raw CSC axes differ')
    if node['indices'].dtype.kind not in 'iu' or node['indices'].dtype.itemsize > 4:
        raise ValueError('This cache contract requires source row indices of at most 32 bits')
    ptr = np.asarray(node['indptr'][:], dtype=np.int64)
    spans = list(column_spans(ptr, max_entries))
    if len(ptr) != shape[1]+1 or len(node['data']) != ptr[-1] or len(node['indices']) != ptr[-1]:
        raise ValueError('CSC arrays cardinality differs')
    contract = dict(source=source_key, rows=axis_hash(row_ids), genes=axis_hash(gene_ids),
                    shape=list(shape), pointer_sha256=hashlib.sha256(ptr.tobytes()).hexdigest(),
                    max_entries=max_entries, spans=[list(x) for x in spans],
                    data_dtype=str(node['data'].dtype), indices_dtype=str(node['indices'].dtype))
    ledger = root/'manifest.json'
    manifest = json.loads(ledger.read_text()) if ledger.exists() else dict(contract=contract, blocks=[])
    if manifest['contract'] != contract:
        raise ValueError('Frozen count source/axis/block contract changed')
    known = {x['file'] for x in manifest['blocks']}
    extras = [p.name for p in root.iterdir() if p.name != 'manifest.json' and p.name not in known]
    if extras:
        raise ValueError('Unaccepted cache files preserved for inspection: '+str(extras[:4]))
    if len(manifest['blocks']) > len(spans):
        raise ValueError('Too many accepted blocks')
    total_bytes = int(ptr[-1])*(node['data'].dtype.itemsize+node['indices'].dtype.itemsize)
    # Dataset handles stay open across all contiguous reads.
    values, indices = node['data'], node['indices']
    completed = 0
    for number, (first, end) in enumerate(spans):
        start, stop = int(ptr[first]), int(ptr[end])
        local_ptr = ptr[first:end+1]-start
        filename = f'block_{number:04d}_{first:05d}_{end:05d}.npz'
        if number < len(manifest['blocks']):
            item = manifest['blocks'][number]
            if (item['first'], item['end'], item['file']) != (first, end, filename):
                raise ValueError('Accepted column coverage changed')
            data, rows, saved_ptr = block_arrays(root, item)
            if not np.array_equal(local_ptr, saved_ptr):
                raise ValueError('Cached CSC pointer differs')
            check_arrays(data, rows, local_ptr, shape[0])
            resumed = True
        else:
            progress('read_raw_counts', block=number, blocks=len(spans),
                     completed_genes=first, total_genes=shape[1],
                     completed_bytes=completed, total_bytes=total_bytes)
            source_data = values[start:stop]
            rows = np.asarray(indices[start:stop], dtype=np.int32)
            check_arrays(source_data, rows, local_ptr, shape[0])
            raw_sha = hashlib.sha256(source_data.tobytes()+rows.tobytes()).hexdigest()
            data = source_data.astype(np.int32)
            path = root/filename
            with (root/(filename+'.partial')).open('xb') as stream:
                np.savez_compressed(stream, data=data, rows=rows, ptr=local_ptr)
                stream.flush()
                import os
                os.fsync(stream.fileno())
            (root/(filename+'.partial')).rename(path)
            item = dict(file=filename, first=first, end=end, nnz=stop-start,
                        source_typed_payload_sha256=raw_sha, bytes=path.stat().st_size,
                        sha256=sha256(path))
            manifest['blocks'].append(item)
            atomic_json(ledger, manifest)
            resumed = False
        completed += (stop-start)*(values.dtype.itemsize+indices.dtype.itemsize)
        progress('raw_counts_checkpointed', block=number+1, blocks=len(spans),
                 completed_genes=end, total_genes=shape[1], completed_bytes=completed,
                 total_bytes=total_bytes, reused=resumed)
    return manifest


def moments(root, manifest, nrows, mt_mask, keep=None):
    """Exact integer row QC and cohort-specific gene sums, using cache only."""
    totals = np.zeros(nrows, dtype=np.int64)
    detected = np.zeros(nrows, dtype=np.int64)
    mt = np.zeros(nrows, dtype=np.int64)
    gene_sums = np.zeros(manifest['contract']['shape'][1], dtype=np.int64)
    for item in manifest['blocks']:
        data, rows, ptr = block_arrays(root, item)
        totals += np.bincount(rows, weights=data, minlength=nrows).astype(np.int64)
        detected += np.bincount(rows, minlength=nrows)
        for j, gene in enumerate(range(item['first'], item['end'])):
            a, b = int(ptr[j]), int(ptr[j+1])
            raw, rr = data[a:b], rows[a:b]
            gene_sums[gene] = raw.sum(dtype=np.int64) if keep is None else raw[keep[rr]].sum(dtype=np.int64)
            if mt_mask[gene]:
                mt[rr] += raw
    return totals, detected, mt, gene_sums


class CountMeans:
    """Exact streamed means exposed to the unchanged original feature rule."""
    def __init__(self, sums, cells):
        self.sums, self.cells = np.asarray(sums, dtype=np.int64), int(cells)
        if self.cells < 1:
            raise ValueError('No analysis cells')

    def mean(self, axis):
        if axis != 0:
            raise ValueError('Only exact gene means available')
        return self.sums.astype(np.float64)/self.cells


def selected_matrix(root, manifest, keep, columns):
    """Build only the original selected raw RNA matrix, preserving both axes."""
    from scipy import sparse
    columns = np.asarray(columns, dtype=np.int64)
    if (np.diff(columns) <= 0).any():
        raise ValueError('Selected source columns must be sorted unique')
    remap = np.full(len(keep), -1, dtype=np.int32)
    remap[keep] = np.arange(int(np.sum(keep)), dtype=np.int32)
    values, indices, pointers = [], [], [0]
    for item in manifest['blocks']:
        selected = columns[(columns >= item['first']) & (columns < item['end'])]
        if not len(selected):
            continue
        data, rows, ptr = block_arrays(root, item)
        for gene in selected:
            j = gene-item['first']
            a, b = int(ptr[j]), int(ptr[j+1])
            rr = remap[rows[a:b]]
            use = rr >= 0
            values.append(data[a:b][use]); indices.append(rr[use])
            pointers.append(pointers[-1]+int(use.sum()))
    return sparse.csc_matrix((np.concatenate(values), np.concatenate(indices),
                              np.asarray(pointers, dtype=np.int64)),
                             shape=(int(np.sum(keep)), len(columns))).tocsr()
