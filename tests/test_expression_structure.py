from pathlib import Path
import sys
import h5py
import numpy as np
from scipy import sparse
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_expression_structure import blocks, normalized_log, centered_pca


def test_streaming_noncontiguous_rows_and_normalization(tmp_path):
    x = np.array([[1, 2, 0], [0, 1, 4], [2, 0, 1], [1, 1, 1]])
    csr = sparse.csr_matrix(x)
    path = tmp_path / 'input.h5'
    with h5py.File(path, 'w') as h:
        g = h.create_group('X')
        g.attrs['encoding-type'] = 'csr_matrix'
        for key in ['data', 'indices', 'indptr']:
            g[key] = getattr(csr, key)
    with h5py.File(path, 'r') as h:
        got = sparse.vstack([part for _, part in blocks(h, np.array([0, 2, 3]), 3, 2)]).toarray()
    np.testing.assert_array_equal(got, x[[0, 2, 3]])
    result = normalized_log(csr, x.sum(axis=1))
    np.testing.assert_allclose(result.toarray(), np.log1p(x / x.sum(axis=1)[:, None] * 10000))
    with pytest.raises(ValueError):
        normalized_log(csr, np.zeros(4))


def test_centered_sparse_pca_matches_dense_svd():
    x = np.random.default_rng(7).poisson(2, size=(20, 8)).astype(float)
    coordinates, loadings, variance = centered_pca(sparse.csr_matrix(x), components=3)
    _, s, vt = np.linalg.svd(x - x.mean(axis=0), full_matrices=False)
    np.testing.assert_allclose(variance, s[:3] ** 2 / 19, rtol=1e-8)
    np.testing.assert_allclose(loadings.T @ loadings, vt[:3].T @ vt[:3], atol=1e-8)
    np.testing.assert_allclose(coordinates.mean(axis=0), 0, atol=1e-8)
