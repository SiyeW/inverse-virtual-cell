import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
import numpy as np
import pytest
from render_target_rna import select_points


def test_plot_keeps_zero_and_low_reference_but_omits_no_reference():
    t=pd.DataFrame(dict(guide=['a','b','c','d'],mapping_status=['exact_catalogue_id']*4,cells=[20,20,20,1],NTC_mean_per_10000=[1,.01,0,1],NTC_detected_fraction=[.5,.001,0,.5],mean_per_10000=[0,.02,0,1],RNA_remaining_ratio=[0,2,np.nan,1]))
    p,good=select_points(t)
    assert p.guide.tolist()==['a','b'] and good.tolist()==[True,False]
    t.loc[0,'RNA_remaining_ratio']=1
    with pytest.raises(ValueError):select_points(t)
