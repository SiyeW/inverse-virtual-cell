import numpy as np
import pytest
from scripts.audit_knockdown_selection import screen


def test_pre_cell_filter_guide_mean_cannot_be_rescued_by_dropping_high_cells():
    r = screen([0, 0, 3, 3], 1, minimum_exclusive=1)
    assert r['ratio_before'] == 1.5
    assert r['ratio_after'] == 0
    assert not r['final_pass']


def test_cell_filter_and_remaining_coverage_match_actual_code_boundaries():
    r = screen([0, 0, 1, 1.5], 1, minimum_exclusive=2)
    assert r['guide_repression_pass']
    assert r['after_cell_gate'] == 3
    assert r['ratio_after'] == pytest.approx(1/3)
    assert r['final_pass']
    assert not screen([0, 0, 1, 1.5], 1, minimum_exclusive=3)['final_pass']


def test_unmeasured_control_stays_undefined_instead_of_a_success():
    assert not screen([0]*30, 0)['final_pass']
    assert not screen([], 1)['reference_defined']
    assert screen([0]*26, .1)['final_pass']  # Arithmetic is not proof of editing.
    with pytest.raises(ValueError):
        screen([0, np.nan], 1)
