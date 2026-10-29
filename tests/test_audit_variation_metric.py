import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_variation_metric import metrics,sample_checks


def test_distances_are_not_correlations_and_shared_center_cancels():
    a=np.array([0.,3.]);b=np.array([3.,0.]);keep=np.ones(2,bool)
    x=metrics(a,b,np.zeros(2),np.zeros(2),keep)
    assert x['mean_absolute_difference']==2 and x['RMS_difference']==2
    assert x['response_pearson']==pytest.approx(-1)
    y=metrics(a,b,np.array([10,20]),np.array([10,20]),keep)
    assert y['mean_absolute_difference']==x['mean_absolute_difference']


def test_common_feature_average_dilutes_a_sparse_difference():
    a=np.zeros(2000);b=a.copy();b[:20]=1
    x=metrics(a,b,a,a,np.ones(2000,bool))
    assert x['mean_absolute_difference']==pytest.approx(.01)
    assert x['RMS_difference']==pytest.approx(.1)
    assert x['response_pearson'] is None


def test_sampling_disjointness_is_checked_instead_of_assumed():
    data=[]
    for n in [5,10]:
        for side in [0,1]:data.append(dict(guide='g',cells=n,repeat=1,side=side,analysis_rows=';'.join(map(str,range(side*10,side*10+n)))))
    p=pd.DataFrame(data);assert sample_checks(p,20)['disjoint_halves_checked']
    p.loc[1,'analysis_rows']='0;11;12;13;14'
    with pytest.raises(ValueError):sample_checks(p,20)
