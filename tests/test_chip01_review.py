from pathlib import Path
import sys,numpy as np,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_chip01_results import descriptive_similarity


def test_descriptive_view_does_not_require_significance_or_assign_groups():
    rng=np.random.default_rng(2);x=rng.normal(size=(5,40));original=x.copy()
    order,corr,columns=descriptive_similarity(x,np.array([10,49,50,80,120]))
    assert set(order)=={2,3,4} and corr.shape==(3,3) and len(columns)==40
    np.testing.assert_allclose(np.diag(corr),1);np.testing.assert_array_equal(x,original)
    order,corr,_=descriptive_similarity(x,np.zeros(5))
    assert len(order)==0 and corr.size==0


def test_descriptive_view_refuses_invalid_profiles():
    with pytest.raises(ValueError):descriptive_similarity(np.array([[np.nan]]),np.array([50]))
