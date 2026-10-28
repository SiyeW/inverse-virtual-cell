from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
import numpy as np
import pytest
from render_target_rna_precision import validate_points,detail_rows


def test_original_point_identity_and_zero_interval_guard():
    old=pd.DataFrame(dict(guide=['a','b'],cells=[10,20],mean_per_10000=[0.,1.],NTC_mean_per_10000=[1.,1.],RNA_remaining_ratio=[0.,1.],NTC_detected_fraction=[.3,.3]))
    new=old.assign(guide_ci_status=['all_zero_no_RNA_mean_interval','estimated'],guide_ci_low=[np.nan,.5],guide_ci_high=[np.nan,1.5],NTC_equal_n_low=[.2,.3])
    validate_points(old,new)
    new.loc[0,'guide_ci_high']=0
    with pytest.raises(ValueError):validate_points(old,new)
    new.loc[0,'guide_ci_high']=np.nan;new.loc[1,'RNA_remaining_ratio']=.5
    with pytest.raises(ValueError):validate_points(old,new)


def test_examples_are_frozen_by_names_and_coverage_not_interval_outcome():
    t=pd.DataFrame(dict(guide=['POU5F1_2','RNF7_1','strongest','RNF7_2'],target=['POU5F1','RNF7','other','RNF7'],cells=[30,30,500,5],NTC_mean_per_10000=[1]*4,RNA_remaining_ratio=[.9,.8,.001,.2]))
    assert detail_rows(t,'kolf').guide.tolist()==['RNF7_1','POU5F1_2']
