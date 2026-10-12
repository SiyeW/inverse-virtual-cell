"""Observable data-contract checks for bounded-memory analysis."""
from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from analyze_guide_capture_qc import measure_cells, join_guides, write_streamed_h5ad, choose_qc


def make_matrix(path):
    # Features x cells, including a custom RNA target, mitochondrion and guide.
    matrix = sparse.csc_matrix(np.array([[4, 2, 5, 1], [1, 0, 2, 1], [0, 3, 0, 1], [7, 0, 9, 2]], dtype=np.int32))
    with h5py.File(path, 'w') as handle:
        g = handle.create_group('matrix')
        for name, value in [('data', matrix.data), ('indices', matrix.indices), ('indptr', matrix.indptr), ('shape', matrix.shape)]:
            g[name] = value
        g['barcodes'] = np.asarray(['a', 'b', 'c', 'd'], dtype='S')
        f = g.create_group('features')
        f['id'] = np.asarray(['ENSG1', 'MT1', 'PB.custom', 'GUIDE'], dtype='S')
        f['name'] = np.asarray(['GENE1', 'MT-CO1', 'PB.custom', 'guide1'], dtype='S')
        f['feature_type'] = np.asarray(['Gene Expression'] * 3 + ['CRISPR Guide Capture'], dtype='S')
    return matrix


def test_full_denominator_and_streamed_counts_round_trip(tmp_path):
    path = tmp_path / 'matrix.h5'
    original = make_matrix(path)
    obs, features, rna = measure_cells(path)
    assert obs['total_counts'].tolist() == [5, 5, 7, 3]
    assert obs['raw_guide_umis'].tolist() == [7, 0, 9, 2]
    calls = pd.DataFrame({'cell_barcode': ['a', 'c', 'd'], 'num_features': [1, 2, 1],
                          'feature_call': ['NTC', 'NTC|PB-guide', 'PB-guide'], 'num_umis': ['7', '8|1', '2']})
    ref = pd.DataFrame({'name': ['NTC', 'PB-guide'], 'target_gene_name': ['Non-Targeting', 'PB.custom']})
    joined = join_guides(obs, calls, ref, 'Non-Targeting')
    assert len(joined) == 4
    assert joined['guide_status'].tolist() == ['single_guide', 'no_confident_call', 'multi_guide', 'single_guide']
    assert joined['perturbed_gene'].tolist() == ['ctrl', '', '', 'PB.custom']
    selected = joined['guide_status'].eq('single_guide').to_numpy()
    out = tmp_path / 'out.h5ad'
    metadata = write_streamed_h5ad(path, joined, features, rna, selected, out, {'test': 'roundtrip'})
    actual = ad.read_h5ad(out)
    np.testing.assert_array_equal(actual.X.toarray(), original[:3, [0, 3]].T.toarray())
    np.testing.assert_array_equal(actual.layers['counts'].toarray(), actual.X.toarray())
    assert actual.obs_names.tolist() == ['a', 'd']
    assert actual.var_names.tolist() == ['ENSG1', 'MT1', 'PB.custom']
    assert metadata['reopened_row_sums_match']


def test_unknown_guide_is_rejected(tmp_path):
    path = tmp_path / 'matrix.h5'
    make_matrix(path)
    obs, _, _ = measure_cells(path)
    calls = pd.DataFrame({'cell_barcode': ['a'], 'num_features': [1], 'feature_call': ['UNKNOWN'], 'num_umis': [2]})
    ref = pd.DataFrame({'name': ['known'], 'target_gene_name': ['GENE1']})
    with pytest.raises(ValueError, match='absent from authoritative'):
        join_guides(obs, calls, ref, 'Non-Targeting')


def test_high_counts_not_silently_called_doublets():
    obs = pd.DataFrame({'total_counts': [1000] * 10 + [100000], 'n_genes_by_counts': [400] * 11,
                        'pct_counts_mt': [1.0] * 11, 'guide_status': ['single_guide'] * 11})
    qc, thresholds = choose_qc(obs, {})
    assert qc.iloc[-1]['high_total_counts']
    assert qc.iloc[-1]['qc_pass_stage1']
    assert qc['doublet_assessment'].eq('pending').all()
