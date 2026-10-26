"""The new adapter must reproduce the adopted method, not just new metrics."""
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import public_count_cache as cache
import validate_original_on_public as challenge
from analyze_guide_capture_qc import choose_qc
from perturbseq_core_io import select_features
from perturbseq_core_effects import effect_profiles, target_support


def source(h, matrix):
    x = sparse.csc_matrix(matrix, dtype=np.float32)
    node = h.create_group('counts')
    node.attrs['encoding-type'] = 'csc_matrix'
    node.attrs['shape'] = x.shape
    for field in ['data', 'indices', 'indptr']:
        node.create_dataset(field, data=getattr(x, field))
    return node


def freeze(node, path, matrix, progress=lambda *_a, **_kw: None):
    return cache.cache_counts(node, path, dict(historical_sha256='known', bytes=123, mtime_ns=5),
                              [f'cell{i}' for i in range(len(matrix))],
                              [f'gene{i}' for i in range(len(matrix[0]))], progress, max_entries=4)


def test_interrupted_read_reuses_only_accepted_blocks_and_restores_without_payload(tmp_path):
    x = np.array([[1, 0, 3, 0, 4], [0, 2, 5, 1, 3], [2, 1, 0, 4, 0]])
    class Interrupted(Exception):
        pass
    def stop(stage, **kw):
        if stage == 'raw_counts_checkpointed':
            raise Interrupted()
    with h5py.File(tmp_path/'source.h5', 'w') as h:
        node = source(h, x)
        with pytest.raises(Interrupted):
            freeze(node, tmp_path/'cache', x, stop)
        first = json.loads((tmp_path/'cache/manifest.json').read_text())['blocks'][0]
        # Poison exactly the accepted source interval. A resume that rereads it
        # will reject fractional counts, rather than accidentally passing.
        cut = int(node['indptr'][first['end']])
        node['data'][:cut] = .25
        manifest = freeze(node, tmp_path/'cache', x)
        # Full cache is independent of the raw payload after publication.
        del node['data']; del node['indices']
        got = cache.selected_matrix(tmp_path/'cache', manifest, np.ones(3, bool), np.arange(5))
        np.testing.assert_array_equal(got.toarray(), x)
        totals, detected, mt, sums = cache.moments(tmp_path/'cache', manifest, 3,
                                                  np.array([True, False, False, False, False]))
        np.testing.assert_array_equal(totals, x.sum(1))
        np.testing.assert_array_equal(detected, (x > 0).sum(1))
        np.testing.assert_array_equal(mt, x[:, 0])
        np.testing.assert_array_equal(sums, x.sum(0))


def test_corruption_and_axis_drift_fail_closed(tmp_path):
    x = np.array([[2, 0, 3], [1, 4, 0]])
    with h5py.File(tmp_path/'source.h5', 'w') as h:
        node = source(h, x)
        manifest = freeze(node, tmp_path/'cache', x)
        with pytest.raises(ValueError, match='contract changed'):
            cache.cache_counts(node, tmp_path/'cache', dict(historical_sha256='changed'),
                               ['cell0', 'cell1'], ['gene0', 'gene1', 'gene2'], lambda *_a, **_kw: None, 4)
        damaged = tmp_path/'cache'/manifest['blocks'][0]['file']
        damaged.write_bytes(damaged.read_bytes()+b'bad')
        with pytest.raises(ValueError, match='content changed'):
            freeze(node, tmp_path/'cache', x)


@pytest.mark.parametrize('bad', ['fraction', 'duplicate', 'negative', 'zero'])
def test_invalid_raw_sparse_values_cannot_enter_adapter(tmp_path, bad):
    x = [[1, 2], [3, 4]]
    with h5py.File(tmp_path/'source.h5', 'w') as h:
        node = source(h, x)
        if bad == 'duplicate':
            node['indices'][1] = node['indices'][0]
        else:
            node['data'][0] = {'fraction': .5, 'negative': -1, 'zero': 0}[bad]
        with pytest.raises(ValueError):
            freeze(node, tmp_path/'cache', x)


def test_adapter_preserves_original_QC_feature_rule_depth_effects_and_target_median(tmp_path):
    rng = np.random.default_rng(123)
    dense = rng.poisson(40, size=(80, 12)).astype(np.int32)
    dense[:, -1] = 0
    dense[:, -2] = 0
    dense[:40, -2] = 1  # Exactly .5 is excluded by the original strict threshold.
    dense[0, :] = 0; dense[0, 0] = 2  # Must be excluded by original RNA QC.
    # ≥200 genes floor would reject a 12-gene toy: extend with expressed genes.
    dense = np.c_[dense, rng.poisson(3, size=(80, 215))].astype(np.int32)
    dense[0] = 0; dense[0, 0] = 2
    with h5py.File(tmp_path/'counts.h5', 'w') as h:
        node = source(h, dense)
        manifest = freeze(node, tmp_path/'cache', dense)
    mt_mask = np.zeros(dense.shape[1], bool)
    totals, detected, mt, _ = cache.moments(tmp_path/'cache', manifest, 80, mt_mask)
    guides = np.repeat(['NTC_1', 'NTC_2', 'GENE_1', 'GENE_2'], 20)
    all_obs = pd.DataFrame(dict(total_counts=totals, n_genes_by_counts=detected,
                                pct_counts_mt=mt/totals*100, guide_id=guides,
                                guide_status='single_guide', source_row=np.arange(80)))
    full, bounds = choose_qc(all_obs, {})
    keep = full.main_analysis_stage1.to_numpy()
    assert not keep[0] and keep[1:].all()
    _, _, _, sums = cache.moments(tmp_path/'cache', manifest, 80, mt_mask, keep)
    var = pd.DataFrame({'gene_name': ['GENE']+[f'rna{i}' for i in range(1, dense.shape[1])]},
                       index=[f'id{i}' for i in range(dense.shape[1])])
    whole = sparse.csr_matrix(dense[keep])
    spec = dict(mode='pooled_raw_mean', raw_mean_gt=.5)
    features = select_features(cache.CountMeans(sums, int(keep.sum())), var, spec)
    np.testing.assert_array_equal(features, select_features(whole, var, spec))
    streamed = cache.selected_matrix(tmp_path/'cache', manifest, keep, features)
    np.testing.assert_array_equal(streamed.toarray(), whole[:, features].toarray())
    obs = full.loc[keep].copy()
    ref = obs.guide_id.str.startswith('NTC').to_numpy()
    scale = obs.loc[ref, 'total_counts'].median()
    norm_streamed = sparse.diags(scale/obs.total_counts.to_numpy(float)) @ streamed.astype(float)
    norm_direct = sparse.diags(scale/obs.total_counts.to_numpy(float)) @ whole[:, features].astype(float)
    overlay = pd.DataFrame(dict(guide=['NTC_1', 'NTC_2', 'GENE_1', 'GENE_2'],
                                 target_gene_name=['Non-Targeting', 'Non-Targeting', 'GENE', 'GENE'],
                                 source_control_type=['non_targeting_control']*2+['targeting']*2))
    arguments = (obs, overlay, ref, var.gene_name.to_numpy()[features], var.gene_name.to_numpy())
    result_a = effect_profiles(norm_streamed.tocsr(), *arguments, seed=20261022)
    result_b = effect_profiles(norm_direct.tocsr(), *arguments, seed=20261022)
    for a, b in zip(result_a, result_b):
        if isinstance(a, pd.DataFrame):
            pd.testing.assert_frame_equal(a, b)
        elif isinstance(a, np.ndarray):
            np.testing.assert_array_equal(a, b)
        else:
            assert a == b
    records, means, refs, effects, first, metrics, selections = result_a
    assert metrics.cells.eq(10).all() and metrics.groupby('guide').size().eq(8).all()
    assert records.reference_overlap_cells.eq(0).all()
    assert np.isnan(effects[records.target.eq('GENE'), 0]).all()
    support, pairs = target_support(norm_streamed.tocsr(), obs, overlay, records, effects, ref)
    assert support.median_pair_r.iloc[0] == pairs.pearson.median()


def test_original_scientific_modules_and_feature_function_are_unchanged():
    proof = challenge.verify_original_functions(ROOT)
    assert len(proof) == 4


def test_background_selection_uses_metadata_and_not_correlation():
    records = pd.DataFrame(dict(guide=['a', 'b', 'c', 'd', 'n1', 'n2'],
                                 cells=[20]*6, target=['T', 'T', 'U', 'V', 'N', 'N'],
                                 source_control_type=['targeting']*4+['non_targeting_control']*2))
    obs = pd.DataFrame(dict(guide_id=np.repeat(records.guide.to_numpy(), 20), batch=['A']*120))
    pairs = pd.DataFrame([dict(target='T', guide_a='a', guide_b='b')])
    effects = np.random.default_rng(9).normal(size=(6, 100))
    first, _ = challenge.diagnostic_backgrounds(obs, records, effects, pairs, 20261022)
    second, _ = challenge.diagnostic_backgrounds(obs, records, -effects[:, ::-1], pairs, 20261022)
    pd.testing.assert_frame_equal(first[['control_a', 'control_b']], second[['control_a', 'control_b']])
    assert first.matched.all()
    assert first.loc[first.comparator.eq('different_target'), 'control_b'].iloc[0] not in ['a', 'b']


@pytest.mark.skipif(sys.platform == 'win32', reason='Actual Linux /proc owner and alarm gate')
def test_entire_original_method_adapter_and_corresponding_plots_on_destination(tmp_path, monkeypatch):
    rng = np.random.default_rng(8)
    values = rng.poisson(5, size=(80, 260)).astype(np.int32)
    # Includes a reproducible same-target response; comparisons remain descriptive.
    values[40:, 10:25] += 20
    names = np.asarray(['GENE']+[f'rna{i}' for i in range(1, 260)], dtype=object)
    guides = np.repeat(['Non-Targeting_1', 'Non-Targeting_2', 'GENE_1', 'GENE_2'], 20)
    path = tmp_path/'public.h5ad'
    with h5py.File(path, 'w') as h:
        layer = h.create_group('layers')
        n = source(layer, values)
        obs = h.create_group('obs'); var = h.create_group('var')
        string = h5py.string_dtype('utf-8')
        for key, arr in dict(_index=[f'c{i}' for i in range(80)], gRNA=guides,
                              gene_target=np.repeat(['NTC', 'NTC', 'GENE', 'GENE'], 20),
                              batch=['ALPHA']*80).items():
            obs.create_dataset(key, data=np.asarray(arr, dtype=object), dtype=string)
        for key, arr in dict(n_gRNA=np.ones(80, int), n_gRNA_UMIs=np.full(80, 5),
                              n_UMI_counts=values.sum(1), n_genes=(values > 0).sum(1)).items():
            obs.create_dataset(key, data=arr)
        var.create_dataset('_index', data=names, dtype=string)
        var.create_dataset('gene_ids', data=np.asarray([f'id{i}' for i in range(260)], object), dtype=string)
        # No dense X at all: the original reader adapter must use counts only.
    cat = tmp_path/'catalogue.csv'
    pd.DataFrame(dict(id=['Non-Targeting-1','Non-Targeting-2','GENE-1','GENE-2'],
                      target_gene_name=['Non-Targeting']*2+['GENE']*2)).to_csv(cat, index=False)
    from perturbseq_core_io import sha256
    config = dict(source=str(path), expected_bytes=path.stat().st_size,
                  expected_mtime_ns=path.stat().st_mtime_ns, historical_sha256=sha256(path),
                  output=str(tmp_path/'result'), allowed_output_parent=str(tmp_path),
                  catalogue=str(cat), catalogue_sha256=sha256(cat),
                  expected_commit='fixture', block_entries=5000, max_seconds=120)
    config_path = tmp_path/'run.json'; config_path.write_text(json.dumps(config))
    monkeypatch.setattr(challenge, 'revision', lambda *_args: 'fixture')
    # The plot hook mutates only presentation globals; isolate it in this test.
    import perturbseq_core_plots as plots
    import run_core_analysis as run
    monkeypatch.setattr(plots, 'save', plots.save)
    monkeypatch.setattr(run, 'save', run.save)
    monkeypatch.setattr(plots, 'NAMES', plots.NAMES.copy())
    monkeypatch.setattr(run, 'NAMES', run.NAMES.copy())
    summary = challenge.run(config_path)
    assert summary['source_cells'] == 80 and summary['response_features'] == 260
    out = tmp_path/'result'
    assert not (out/'owner.lock').exists()
    assert json.loads((out/'status.json').read_text())['stage'] == 'complete_pending_method_interpretation'
    for name in ['01_technical_qc','02_four_type_coverage','03_four_type_cell_quality',
                 '05_matched_rna_effects','07_target_support']:
        assert (out/'figures'/f'{name}.png').stat().st_size > 1000
    target = pd.read_csv(out/'tables/within_target_pairs.csv')
    assert target.pearson.iloc[0] > .8
