import sys
from pathlib import Path
import numpy as np,pandas as pd,pytest
from scipy import sparse
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from audit_scrublet_workflow import prepare_native_inputs


def fixture_counts():
    rng=np.random.default_rng(4)
    means=rng.gamma(1,1,size=(240,800))
    means[:120,:100]*=5;means[120:,100:200]*=5
    counts=sparse.csr_matrix(rng.poisson(means))
    obs=pd.DataFrame(index=[f'cell{i}' for i in range(240)])
    var=pd.DataFrame(index=[f'gene{i}' for i in range(800)])
    return counts,obs,var


def test_native_preparation_preserves_raw_and_ignores_guide_labels():
    counts,obs,var=fixture_counts();before=counts.copy()
    obs['guide_id']=['a']*120+['b']*120
    x,sim,info=prepare_native_inputs(counts,obs,var)
    swapped=obs.copy();swapped['guide_id']=obs.guide_id.iloc[::-1].to_numpy()
    y,_,other=prepare_native_inputs(counts,swapped,var)
    assert (counts!=before).nnz==0 and x.obs_names.equals(obs.index)
    assert x.var_names.equals(y.var_names) and (x.X!=y.X).nnz==0
    assert info==other and info['hvg_selection_passes']==1
    assert 30<x.n_vars<800 and sim.n_obs==480


def test_scores_match_native_wrapper_on_same_raw_input():
    import anndata as ad,scanpy as sc
    from threadpoolctl import threadpool_limits
    counts,obs,var=fixture_counts()
    raw=ad.AnnData(counts.copy(),obs=obs.copy(),var=var.copy())
    prepared,sim,_=prepare_native_inputs(counts,obs,var,seed=0)
    settings=dict(expected_doublet_rate=.05,sim_doublet_ratio=2.,mean_center=True,n_prin_comps=10,use_approx_neighbors=False,random_state=0)
    with threadpool_limits(limits=1):
        sc.pp.scrublet(raw,**settings)
        sc.pp.scrublet(prepared,adata_sim=sim,**settings)
    np.testing.assert_allclose(raw.obs.doublet_score,prepared.obs.doublet_score,rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(raw.uns['scrublet']['doublet_scores_sim'],prepared.uns['scrublet']['doublet_scores_sim'],rtol=1e-12,atol=1e-12)


def test_invalid_counts_are_rejected():
    counts,obs,var=fixture_counts();counts.data[0]=-1
    with pytest.raises(ValueError,match='nonnegative integer'):
        prepare_native_inputs(counts,obs,var)
