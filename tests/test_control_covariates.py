from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from review_control_covariates import residualize


def test_residuals_orthogonal_and_constant_covariate_safe():
    rng=np.random.default_rng(3)
    z=np.column_stack([rng.normal(size=50),np.ones(50)])
    y=z[:,:1]*np.array([[2,3]])+rng.normal(size=(50,2))
    residual,rank=residualize(y,z)
    assert rank==2
    np.testing.assert_allclose(residual.mean(0),0,atol=1e-12)
    np.testing.assert_allclose(z[:,0]@residual,0,atol=1e-12)
    with pytest.raises(ValueError):
        residualize(y,np.array([[np.nan]]))
