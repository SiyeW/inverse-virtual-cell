import sys
from pathlib import Path
import numpy as np
import pytest
from scipy import sparse
pytest.importorskip('scanpy')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from recover_frozen_rna_counts import recover_counts,_forward


def test_independent_full_gene_normalization_then_panel_selection():
    import anndata as ad
    import scanpy as sc
    rng=np.random.default_rng(51);raw=rng.poisson(3,size=(32,19)).astype('int32')
    raw[0,2]=800;totals=raw.sum(axis=1)
    a=ad.AnnData(sparse.csr_matrix(raw.astype('float32')))
    sc.pp.normalize_total(a,target_sum=10000);sc.pp.log1p(a)
    selected=np.array([1,2,6,11,18]);panel=a.X[:,selected].tocsr()
    recovered,proof=recover_counts(panel,totals,chunk_entries=17)
    assert np.array_equal(recovered.toarray(),raw[:,selected]) and proof['adjacent_integer_uniqueness']


def test_known_panel_counts_recover_exactly_with_omitted_rna_totals():
    counts=sparse.csr_matrix(np.array([[0,1,7],[5000,0,31],[12,21,2]],dtype='int32'))
    totals=np.array([200,17000,700])
    accepted=_forward(counts,totals,10000)
    recovered,proof=recover_counts(accepted,totals,chunk_entries=3)
    assert np.array_equal(recovered.toarray(),counts.toarray())
    assert proof['roundtrip_bit_exact'] and proof['adjacent_integer_uniqueness']


def test_noninvertible_edited_cache_fails_instead_of_rounding_silently():
    counts=sparse.csr_matrix([[2,5],[7,4]],dtype='int32');totals=np.array([100,200])
    accepted=_forward(counts,totals,10000);accepted.data[0]+=.013
    with pytest.raises(ValueError):recover_counts(accepted,totals)


def test_invalid_total_axis_and_precision_are_rejected():
    counts=sparse.csr_matrix([[2,5],[7,4]],dtype='int32');accepted=_forward(counts,np.array([100,200]),10000)
    with pytest.raises(ValueError):recover_counts(accepted,np.array([100]))
    with pytest.raises(ValueError):recover_counts(accepted.astype('float64'),np.array([100,200]))
