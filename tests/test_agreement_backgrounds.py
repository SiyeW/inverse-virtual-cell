from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import pytest
from refine_agreement_backgrounds import matched_coverage,pooled_control_profiles

def test_controls_exclude_whole_guide_across_batches():
    controls=pd.DataFrame(dict(guide=['A','A','B','B'],cells=[1,1,1,1]))
    source=np.array([[3.],[5.],[10.],[12.]])
    old=pooled_control_profiles(controls,source,'old')
    assert old[0,0]==pytest.approx(np.log2(4/12))
    assert old[1,0]==pytest.approx(np.log2(6/12))
    logged=np.log2(source+1)
    new=pooled_control_profiles(controls,logged,'new')
    assert new[0,0]==pytest.approx(np.log2(4)-(np.log2(11)+np.log2(13))/2)

def test_count_matching_is_fixed_before_effects():
    plan=pd.DataFrame(dict(cells_a=[10,10],cells_b=[20,20],NTC_cells_a=[12,100],NTC_cells_b=[25,20],other_cells_a=[9,10],other_cells_b=[21,20]))
    assert list(matched_coverage(plan))==[True,False]

def test_no_independent_reference_fails():
    with pytest.raises(ValueError):pooled_control_profiles(pd.DataFrame(dict(guide=['A','A'],cells=[1,1])),np.ones((2,3)),'old')
