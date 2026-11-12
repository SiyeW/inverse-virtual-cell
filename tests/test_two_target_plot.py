"""The RNA pilot plot retains low-count guides and missing-score states."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
from plot_two_target_rna import draw


def test_plot_keeps_all_finite_records_and_unscored_guide(tmp_path):
    rng = np.random.default_rng(41)
    rows = []
    names = ['NTC', 'No sgRNA detected', 'Intergenic'] + [f'Gene-{i}' for i in range(1, 6)]
    for name, number in zip(names, [50, 30, 10, 18, 24, 12, 16, 20]):
        values = rng.normal(size=number)
        if name == 'Gene-5':
            values[:] = np.nan
        frame = pd.DataFrame({'score': values, 'indel_support': np.zeros(number, dtype=bool)})
        if name == 'Gene-1':
            frame.loc[0, 'indel_support'] = True
        rows.append((name, frame))
    counts = pd.DataFrame({'guide_id': names[3:], 'DE_RNAs': [0, 2, 3, 0, 1],
                           'assessment': ['evaluated'] * 5}).set_index('guide_id')
    result = draw('Gene', rows, counts, (-4, 4), 4, tmp_path / 'fixture_only.png')
    assert result['finite_score_records'] == sum(len(frame) for _, frame in rows[:-1])
    assert result['plotted_indel_records'] == 1
    assert result['rows'][-1] == {'label': 'Gene-5', 'records': 20, 'scored': 0}
    assert len(result['rows']) == 8


def test_absent_model_never_plots_fallback_scores(tmp_path):
    rows = [(name, pd.DataFrame({'score': [np.nan] * 3, 'indel_support': [True, False, False]}))
            for name in ['NTC', 'No sgRNA detected', 'Intergenic'] + [f'Gene-{i}' for i in range(1, 6)]]
    counts = pd.DataFrame({'guide_id': [name for name, _ in rows[3:]], 'DE_RNAs': [0, 0, 1, 0, 2],
                          'assessment': ['evaluated'] * 5}).set_index('guide_id')
    result = draw('Gene', rows, counts, (-3, 3), 3, tmp_path / 'unscored_fixture.png', feature_count=1)
    assert result['score_model'] == 'not_estimable'
    assert result['finite_score_records'] == 0 and result['plotted_indel_records'] == 0
    assert result['indel_support_records_shown'] == 5
    assert len(result['rows']) == 8
