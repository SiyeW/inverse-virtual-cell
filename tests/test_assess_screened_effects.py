import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_screened_effects import revised_mean, pair_plan, reference_row


def test_subtraction_is_exactly_the_kept_population_mean():
    x=np.array([[0,2,4],[2,1,0],[5,2,2],[1,0,3]],float)
    actual=revised_mean(x.mean(0),4,x[[1,3]].sum(0),2)
    np.testing.assert_allclose(actual,x[[0,2]].mean(0),rtol=0,atol=1e-14)
    with pytest.raises(ValueError):revised_mean([0,1],3,[1,0],2)


def test_reference_prefix_reads_the_named_row_without_requiring_whole_array(tmp_path):
    p=tmp_path/'profiles.npz';a=np.arange(30,dtype=float).reshape(5,6)
    np.savez_compressed(p,reference_means=a)
    x,proof=reference_row(p,1)
    np.testing.assert_array_equal(x,a[1]);assert not proof['whole_member_crc_checked']
    assert proof['prefix_payload_bytes']==2*6*8


def test_background_plan_is_coverage_based_disjoint_and_reproducible():
    d=pd.DataFrame(dict(guide=['a1','a2','b1','b2','c1'],target=['a','a','b','b','c'],
                        cells=[30,60,31,59,40],source_control_type=['targeting']*5))
    p=pair_plan(d,set(d.guide))
    pd.testing.assert_frame_equal(p,pair_plan(d.sample(frac=1,random_state=3),set(d.guide)))
    for r in p.itertuples():
        assert not {r.guide_a,r.guide_b}&{r.other_a,r.other_b}
        targets=d.set_index('guide').target
        assert targets[r.other_a]!=targets[r.other_b]
