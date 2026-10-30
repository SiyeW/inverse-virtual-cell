"""Recover panel counts only when the recorded float32 transform is invertible.

Authoritative full-RNA totals are required. Every recovered entry must recreate
the stored normalized/log value bit-for-bit, and adjacent valid integers must
produce different values. This is not treating normalized values as counts.
"""
import numpy as np
from scipy import sparse


def _forward(counts, totals, target_sum):
    import anndata as ad
    import scanpy as sc
    values=counts.astype(np.float32)
    remaining=np.asarray(totals,dtype=np.int64)-np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    if np.any(remaining<0):raise ValueError('Panel counts exceed authoritative full-RNA totals')
    augmented=sparse.hstack([values,sparse.csr_matrix(remaining[:,None].astype(np.float32))],format='csr')
    a=ad.AnnData(augmented)
    sc.pp.normalize_total(a,target_sum=target_sum)
    sc.pp.log1p(a)
    result=a.X[:,:counts.shape[1]].tocsr();result.eliminate_zeros();result.sort_indices()
    return result


def recover_counts(log_values,full_rna_totals,target_sum=10000,chunk_entries=131072):
    x=sparse.csr_matrix(log_values,copy=True);x.eliminate_zeros();x.sort_indices()
    if x.dtype!=np.float32 or not np.isfinite(x.data).all() or np.any(x.data<0):raise ValueError('Expected finite nonnegative float32 log checkpoint')
    totals=np.asarray(full_rna_totals)
    if totals.shape!=(x.shape[0],) or np.any(totals<=0) or np.any(totals>2**20) or np.any(totals!=np.rint(totals)):
        raise ValueError('Authoritative totals outside exact-integer verification range')
    rows=np.repeat(np.arange(x.shape[0]),np.diff(x.indptr))
    inverse=np.expm1(x.data.astype(np.float64))*totals[rows]/target_sum
    candidate=np.rint(inverse)
    if np.any(candidate<1) or np.any(candidate>totals[rows]):raise ValueError('Invalid recovered count candidate')
    counts=sparse.csr_matrix((candidate.astype(np.int32),x.indices.copy(),x.indptr.copy()),shape=x.shape)
    forward=_forward(counts,totals,target_sum)
    if not np.array_equal(forward.indptr,x.indptr) or not np.array_equal(forward.indices,x.indices) or not np.array_equal(forward.data,x.data):
        raise ValueError('Count recovery does not exactly reproduce the accepted transform')
    # The transform is monotone. Strictly different adjacent valid integers,
    # evaluated by the same public normalization, establish local uniqueness.
    for start in range(0,x.nnz,chunk_entries):
        stop=min(start+chunk_entries,x.nnz);c=candidate[start:stop].astype(np.int64);t=totals[rows[start:stop]].astype(np.int64)
        lower=c-1;upper=np.minimum(c+1,t)
        probe_counts=np.concatenate([lower,upper]);probe_totals=np.concatenate([t,t])
        probe=sparse.csr_matrix(probe_counts[:,None])
        values=_forward(probe,probe_totals,target_sum).toarray().ravel()
        k=len(c);observed=x.data[start:stop]
        if not np.all(values[:k]<observed) or not np.all(values[k:][c<t]>observed[c<t]):
            raise ValueError('Adjacent integer ambiguity; original count subset is required')
    return counts,dict(entries=x.nnz,shape=list(x.shape),roundtrip_bit_exact=True,adjacent_integer_uniqueness=True,
                       max_inverse_rounding_residual=float(np.max(np.abs(inverse-candidate),initial=0)),
                       authoritative_totals_required=True,source='recovered subset from accepted invertible normalized/log checkpoint; not a fresh original-container hash')
