from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]/'scripts'))
from export_original_plot_data import curves
from render_original_method_comparison import render
import json


def inputs():
    records = pd.DataFrame(dict(guide=['N1','N2','T1'],
        source_control_type=['non_targeting_control','non_targeting_control','targeting']))
    # N1 has one finite RNA, N2 has two; each GUIDE must have equal total mass.
    matched = np.array([[1.,np.nan],[0.,0.],[.2,.4]])
    metrics = pd.DataFrame([dict(guide=g, repeat=r, cells=10,
                                 rms_log2fc={'N1':.1,'N2':.3,'T1':.2}[g])
                             for g in records.guide for r in range(1,9)])
    support = pd.DataFrame(dict(target=['T'], source_control_type=['targeting'],
        eligible_guides=[2], assessment=['descriptive_assessable'],
        largest_guide_fraction=[.7], median_pair_r=[.25], pairs=[1]))
    return records, matched, metrics, support


def test_original_guide_weighting_and_exact_ecdf_target_aggregation():
    payload = curves(*inputs())
    ntc = payload['roles']['non_targeting_control']
    edges = np.asarray(payload['density_edges'])
    density = np.asarray(ntc['density'])
    np.testing.assert_allclose(density@np.diff(edges), 1.)
    assert np.isclose((density*np.diff(edges))[edges[:-1]>=.9].sum(), .5)
    assert ntc['rms_x'] == [.1,.1,.3] and ntc['rms_y'] == [0.,.5,1.]
    assert payload['agreement']['targeting']['median_pair_r'] == [.25]
    assert payload['target_coverage'] == [0,0,1,0]
    assert payload['roles']['targeting']['rms_singleton']


@pytest.mark.parametrize('bad', ['cells','repeat','guide_axis'])
def test_changed_sampling_contract_or_guide_axis_is_rejected(bad):
    records, matched, metrics, support = inputs()
    if bad == 'cells': metrics.loc[0,'cells'] = 25
    elif bad == 'repeat': metrics.loc[0,'repeat'] = 2
    else: records.loc[1,'guide'] = 'N1'
    with pytest.raises(ValueError): curves(records, matched, metrics, support)


def test_comparison_requires_public_original_result_and_preserves_facts(tmp_path):
    payload = curves(*inputs())
    payload.update(dataset_id='lab', label='Lab fixture', analysis_cells=30,
                   matched_cells=10, matched_repeats=8,
                   agreement_unit='one_point_per_target_median_of_guide_pairs')
    a=tmp_path/'a.json';a.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='KOLF'):render([a],tmp_path/'missing')
    payload.update(dataset_id='kolf_original',label='Public fixture')
    b=tmp_path/'b.json';b.write_text(json.dumps(payload))
    result=render([a,b],tmp_path/'plots')
    table=pd.read_csv(tmp_path/'plots/tables/original_comparison_facts.csv')
    assert table.median_target_r.eq(.25).all()
    assert table.NTC_RMS.eq(.2).all()
    assert len(list((tmp_path/'plots/figures').glob('*.png')))==2
    payload['matched_cells']=25;b.write_text(json.dumps(payload))
    with pytest.raises(ValueError,match='Changed method'):render([a,b],tmp_path/'wrong_method')
