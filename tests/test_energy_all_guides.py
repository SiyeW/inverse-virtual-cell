import sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from energy_all_frozen_guides import estimates

def test_all_records_and_whole_own_NTC_guide_exclusion():
    rng=np.random.default_rng(75);pc=rng.normal(size=(57,5))
    guide=np.repeat(['nt1','nt2','a','b'],[10,14,13,20]);role=np.where(np.char.startswith(guide,'nt'),'non_targeting_control','targeting')
    got=estimates(pc,guide,role).set_index('guide')
    for g in np.unique(guide):
        q=guide==g;r=(role=='non_targeting_control') & (guide!=g)
        x=pc[q];y=pc[r];n=len(x);m=len(y)
        want=2*cdist(x,y).mean()-cdist(x,x).sum()/(n*(n-1))-cdist(y,y).sum()/(m*(m-1))
        assert abs(got.loc[g,'energy_unbiased']-want)<1e-12
        assert got.loc[g,'records']==n and got.loc[g,'reference_records']==m
