from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from perturbation_statistics import energy_v_statistic, energy_permutation, differential_expression_block, bh_adjust


def test_energy_value_and_symmetry():
    x, y = np.zeros((5, 1)), np.ones((8, 1))
    assert energy_v_statistic(x, y) == 2
    assert energy_v_statistic(y, x) == 2
    assert energy_v_statistic(x, x) == 0


def test_permutation_null_and_separation_reproducible():
    x = np.zeros((20, 2))
    assert energy_permutation(x, x, permutations=99) == (0, 1)
    assert energy_permutation(x, x + 10, permutations=99)[1] == .01
    assert energy_permutation(x, x + 10, permutations=99) == energy_permutation(x, x + 10, permutations=99)


def test_memory_and_invalid_inputs_rejected():
    with pytest.raises(MemoryError):
        energy_permutation(np.ones((20, 2)), np.ones((20, 2)), max_pairwise_mib=.001)
    with pytest.raises(ValueError):
        energy_v_statistic(np.array([[np.nan]]), np.zeros((1, 1)))
    with pytest.raises(ValueError):
        energy_v_statistic(np.zeros((0, 1)), np.zeros((1, 1)))


def test_de_block_ties_effects_and_global_bh():
    x = np.column_stack([np.zeros(20), np.ones(20) * 3, np.arange(20)])
    y = np.column_stack([np.zeros(30), np.zeros(30), np.arange(30)])
    r = differential_expression_block(x, y)
    assert r['pvalue'][0] == 1
    assert r['rank_biserial'][0] == 0
    assert r['rank_biserial'][1] == 1
    assert r['log2_fold_change_pseudocount1'][1] == 2
    np.testing.assert_allclose(r['pvalue'], [mannwhitneyu(x[:,i],y[:,i],method='asymptotic').pvalue for i in range(3)])
    np.testing.assert_allclose(bh_adjust([.01,.04,.02]), [.03,.04,.03])
    with pytest.raises(ValueError):
        differential_expression_block(x[:5], y)
