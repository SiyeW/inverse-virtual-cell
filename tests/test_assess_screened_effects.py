import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import h5py
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assess_screened_effects import revised_mean, pair_plan, reference_row, read_removed, source_identity


def test_subtraction_is_exactly_the_kept_population_mean():
    x=np.array([[0,2,4],[2,1,0],[5,2,2],[1,0,3]],float)
    actual=revised_mean(x.mean(0),4,x[[1,3]].sum(0),2)
    np.testing.assert_allclose(actual,x[[0,2]].mean(0),rtol=0,atol=1e-14)
    with pytest.raises(ValueError):revised_mean([0,1],3,[1,0],2)


def test_reference_prefix_reads_the_named_row_without_requiring_whole_array(tmp_path):
    p=tmp_path/'profiles.npz';a=np.arange(30,dtype=float).reshape(5,6)
    np.savez_compressed(p,reference_means=a)
    x,proof=reference_row(p,1)
    np.testing.assert_array_equal(x,a[1]);assert not proof['whole_member_crc_checked']
    assert proof['prefix_payload_bytes']==2*6*8


def test_background_plan_is_coverage_based_disjoint_and_reproducible():
    d=pd.DataFrame(dict(guide=['a1','a2','b1','b2','c1'],target=['a','a','b','b','c'],
                        cells=[30,60,31,59,40],source_control_type=['targeting']*5))
    p=pair_plan(d,set(d.guide))
    pd.testing.assert_frame_equal(p,pair_plan(d.sample(frac=1,random_state=3),set(d.guide)))
    for r in p.itertuples():
        assert not {r.guide_a,r.guide_b}&{r.other_a,r.other_b}
        targets=d.set_index('guide').target
        assert targets[r.other_a]!=targets[r.other_b]


def test_h5ad_categorical_axes_and_checkpoint_resume_skip_data(tmp_path,monkeypatch):
    source=tmp_path/'cells.h5ad'
    with h5py.File(source,'w') as f:
        for name,values,index in [('obs',['bc1','bc2'],'cell_barcode'),('var',['rna1','rna2'],'feature_id')]:
            group=f.create_group(name);group.attrs['_index']=index
            cat=group.create_group(index)
            cat.create_dataset('categories',data=np.asarray(values,dtype=h5py.string_dtype()))
            cat.create_dataset('codes',data=[0,1])
        node=f.create_group('layers/counts');node.attrs['encoding-type']='csr_matrix'
        node.create_dataset('data',data=[3,1,2]);node.create_dataset('indices',data=[0,1,1]);node.create_dataset('indptr',data=[0,2,3])
    config=dict(path=str(source),identity=source_identity(source),format='csr_h5ad_counts')
    cells=pd.DataFrame(dict(cell_barcode=['bc2','bc1'],total_counts=[2,4]))
    genes=pd.DataFrame(dict(gene_id=['rna2','rna1']))
    request=[dict(row=1,target_column=1,target_count=3)]
    rows,proof=read_removed(config,cells,genes,request,tmp_path/'cache',lambda *a,**k:None)
    np.testing.assert_array_equal(rows[1].toarray(),[[1,3]])
    old=h5py.Dataset.__getitem__
    def refuse_data(self,key):
        if self.name in ['/layers/counts/data','/layers/counts/indices']:raise AssertionError('accepted payload reread')
        return old(self,key)
    monkeypatch.setattr(h5py.Dataset,'__getitem__',refuse_data)
    resumed,proof=read_removed(config,cells,genes,request,tmp_path/'cache',lambda *a,**k:None)
    np.testing.assert_array_equal(resumed[1].toarray(),[[1,3]])
    assert proof['logical_array_bytes']==0
