from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
import pytest
from render_guide_variation import validate_points,DATASETS,KINDS

def fixture(root):
    for name in DATASETS:
        (root/name).mkdir()
        pd.DataFrame([dict(pair_id=i,comparison=k,cells_per_group=10,mean_absolute_difference=.1) for i in [0,1] for k in KINDS]).to_csv(root/name/'variation_points.csv',index=False)

def test_fixed_four_arm_cohort(tmp_path):
    fixture(tmp_path);assert len(validate_points(tmp_path,10))==3

def test_missing_control_or_negative_distance_fails(tmp_path):
    fixture(tmp_path);p=tmp_path/'company/variation_points.csv';t=pd.read_csv(p);t=t.iloc[:-1];t.to_csv(p,index=False)
    with pytest.raises(ValueError):validate_points(tmp_path,10)
