from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import pandas as pd
import pytest
from render_agreement_overlay import load_tables,DATASETS

def tables(root):
    for name in DATASETS:
        for method in ['old','new']:
            pd.DataFrame([dict(target='GENE',largest_guide_fraction=.5,median_pair_r=.2 if method=='old' else .1)]).to_csv(root/f'{name}_{method}_points.csv',index=False)

def test_old_new_same_cells_and_targets(tmp_path):
    tables(tmp_path);result=load_tables(tmp_path)
    assert len(result)==6 and result['kolf','new'].median_pair_r.iloc[0]==.1

@pytest.mark.parametrize('bad',['target','fraction','singleton'])
def test_changed_cohort_or_singleton_is_rejected(tmp_path,bad):
    tables(tmp_path);p=tmp_path/'company_new_points.csv';t=pd.read_csv(p)
    if bad=='target':t.loc[0,'target']='OTHER'
    elif bad=='fraction':t.loc[0,'largest_guide_fraction']=.6
    else:t.loc[0,'largest_guide_fraction']=1.
    t.to_csv(p,index=False)
    with pytest.raises(ValueError):load_tables(tmp_path)
