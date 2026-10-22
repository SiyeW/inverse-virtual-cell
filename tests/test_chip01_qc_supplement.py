from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parents[1]/'scripts'))
from chip01_qc_supplement import (distribution_summary, validate_guide_coverage,
    plan_null_splits, library_normalize, masked_response_profiles,
    profile_correlations, guide_agreement, null_test_summary)


def test_zero_guides_retained_no_pseudocount_in_skew():
    actual = distribution_summary([0, 0, 1, 2, 8], [1, 5])
    assert actual['n'] == 5 and actual['zero_n'] == 2
    assert actual['p90_over_p10'] is None and actual['n_below_1'] == 2
    assert distribution_summary([0, 0])['gini'] is None
    with pytest.raises(ValueError): distribution_summary([1, np.nan])


def test_coverage_requires_full_source_assignment_identity():
    obs = pd.DataFrame(dict(guide_status=['single_guide']*3, guide_id=['nt', 'x', 'x'],
                            main_analysis_stage1=[True, True, False]))
    ref = pd.DataFrame(dict(name=['nt', 'x', 'absent'], target_gene_name=['Non-Targeting', 'G', 'G'],
                           single_guide_before_qc=[1, 2, 0], single_guide_after_stage1_qc=[1, 1, 0]))
    result = validate_guide_coverage(obs, ref)
    assert len(result) == 3 and list(result.guide_type) == ['NTC', 'targeting', 'targeting']
    bad = ref.copy(); bad.loc[0, 'single_guide_before_qc'] = 2
    with pytest.raises(ValueError): validate_guide_coverage(obs, bad)


def test_null_splits_have_no_internal_overlap_and_reproduce():
    first = plan_null_splits(30, [1, 5], 20, 3, 12)
    second = plan_null_splits(30, [1, 5], 20, 3, 12)
    assert len(first) == 6
    for a, b in zip(first, second):
        assert not np.intersect1d(a['treatment'], a['control']).size
        assert len(a['control']) == 20 and len(a['treatment']) == a['cells']
        np.testing.assert_array_equal(a['treatment'], b['treatment'])
        np.testing.assert_array_equal(a['control'], b['control'])
    with pytest.raises(ValueError): plan_null_splits(20, [5], 20, 2, 0)


def test_normalization_preserves_raw_counts_and_rejects_invalid_totals():
    original = sparse.csr_matrix([[1, 3], [4, 4]])
    before = original.copy()
    result = library_normalize(original, [4, 8], 10)
    np.testing.assert_allclose(result.toarray().sum(1), [10, 10])
    assert (original != before).nnz == 0
    with pytest.raises(ValueError): library_normalize(original, [0, 8], 10)


def test_guide_profiles_mask_only_unambiguous_direct_rna_not_guessed_host():
    x = sparse.csr_matrix([[4, 1, 2], [4, 1, 2], [1, 4, 3], [1, 4, 3], [2, 3, 5]])
    rec, profiles = masked_response_profiles(x, ['a', 'a', 'b', 'b', 'low'],
        ['G', 'G', 'PB_ORF', 'PB_ORF', 'G'], np.ones(3), np.ones(3), ['G', 'H', 'I'], 2)
    assert list(rec.guide) == ['a', 'b']
    assert profiles[0, 0] == 0 and profiles[1, 1] == 3
    assert rec.direct_target_feature_masked.tolist() == [True, False]
    with pytest.raises(ValueError):
        masked_response_profiles(x, ['a']*5, ['G', 'H', 'G', 'G', 'G'], np.ones(3), np.ones(3), ['G', 'H', 'I'], 2)


def test_constant_profiles_explicitly_unassessable():
    corr, valid = profile_correlations([[1, 1, 1], [0, 1, 3], [0, 2, 6]])
    assert valid.tolist() == [False, True, True]
    assert np.isnan(corr[0]).all()
    np.testing.assert_allclose(corr[1, 2], 1)


def test_pair_types_do_not_mix_ntc_with_targeting_or_fake_replicates():
    records = pd.DataFrame(dict(guide=list('abcdef'), target=['G', 'G', 'H', 'H', 'ctrl', 'ctrl'],
                                cells=[10, 12, 20, 30, 15, 18],
                                guide_type=['targeting']*4+['NTC']*2))
    profiles = np.random.default_rng(9).normal(size=(6, 20))
    pairs, _ = guide_agreement(records, profiles, max_between=2, seed=1)
    assert pairs.pair_type.value_counts().to_dict() == {'same-target': 2, 'different-target': 2, 'NTC-NTC': 1}


def test_bh_is_once_over_complete_gene_family_not_per_block():
    result = null_test_summary([.001, .02, .04, .9])
    assert result['tested_features'] == 4 and result['bh_q_lt_005'] == 2
    with pytest.raises(ValueError): null_test_summary([np.nan, .1])


def test_real_io_null_smoke_retains_nan_untested_and_full_family(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    from run_chip01_qc_supplement import negative_control_checks
    rng = np.random.default_rng(35)
    counts = sparse.csr_matrix(rng.poisson(.7, (350, 8)))
    # All source totals are positive and can also contain untested features.
    totals = np.asarray(counts.sum(1)).ravel()+10
    obs = pd.DataFrame(dict(is_control=[True]*300+[False]*50, total_counts=totals),
                       index=pd.Index([f'cell{i}' for i in range(350)], dtype='string'))
    genes = pd.DataFrame(dict(relative_z_defined=[True]*8), index=[f'gene{i}' for i in range(8)])
    core = np.arange(350) < 200
    calibration = tmp_path/'calibration.csv'
    pd.DataFrame(dict(statistic=np.linspace(-10, 100, 100),
                      pvalue=np.exp(-np.linspace(0, 15, 100)))).to_csv(calibration,index=False)
    output=tmp_path/'result';output.mkdir()
    (output/'tables').mkdir();(output/'figures').mkdir()
    runtime=tmp_path/'runtime';runtime.mkdir()
    states=[]
    config=dict(expected_tested_genes=8,calibration=str(calibration),seed=4,
                null_gene_block_size=3,null_sizes=[1,5,25],null_repeats=2,null_reference_cells=200)
    result=negative_control_checks(counts,obs,genes,core,config,output,runtime,
                                  lambda stage,**values:states.append((stage,values)),True)
    assert result['pseudo_comparisons']==4 and result['tested_genes']==8
    saved=np.load(output/'tables/null_pvalues.npz')
    assert saved['ad_pvalues'].shape==(4,8)
    assert np.isnan(saved['mw_pvalues'][:2]).all()
    assert np.isfinite(saved['mw_pvalues'][2:]).all()
    rows=pd.read_csv(output/'tables/null_comparison_summary.csv')
    assert len(rows)==6 and rows.tested_features.eq(8).all()
    assert states[-1][1]['completed_genes']==8
