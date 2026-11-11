"""Regressions for whole-gene count transport and source identifier preservation."""
import gzip
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from prepare_two_target_rna import write_sparse_counts


def test_sparse_vectors_are_exact_R_CSC_storage(tmp_path):
    cells_by_genes = sparse.csr_matrix(np.array([
        [0, 2, 7, 0], [3, 0, 0, 0], [0, 4, 0, 9],
    ], dtype=np.int32))
    files = write_sparse_counts(cells_by_genes, tmp_path)
    vectors = []
    for name in files:
        with gzip.open(tmp_path / name, 'rb') as stream:
            vectors.append(np.frombuffer(stream.read(), dtype='<i4'))
    genes_by_cells = sparse.csc_matrix(tuple(vectors), shape=(4, 3))
    np.testing.assert_array_equal(genes_by_cells.toarray(), cells_by_genes.T.toarray())
    np.testing.assert_array_equal(np.asarray(genes_by_cells.sum(axis=0)).ravel(), [9, 3, 13])
    assert not list(tmp_path.glob('*.partial'))
