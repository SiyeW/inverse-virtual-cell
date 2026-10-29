import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_official_mixscape import checkpoint_strings, model_labels, summarize_guides
import run_official_mixscape as workflow


def test_model_groups_keep_intergenic_separate():
    cells = pd.DataFrame(dict(target=['A', 'A', '', '', ''],
                              guide_id=['A-1', 'A-2', 'NT-1', 'IG-1', 'IG-2'],
                              role=['targeting', 'targeting', 'non_targeting_control',
                                    'intergenic_control', 'intergenic_control']))
    assert list(model_labels(cells)) == ['A', 'A', 'NT', 'INTERGENIC:IG-1', 'INTERGENIC:IG-2']


def test_missing_guides_and_NTC_are_not_zero_efficiency():
    catalogue = pd.DataFrame(dict(guide=['A-1', 'A-2', 'B-1', 'NT-1'],
                                  target_gene_name=['A', 'A', 'B', ''],
                                  source_control_type=['targeting'] * 3 + ['non_targeting_control']))
    cells = pd.DataFrame(dict(cell_barcode=['c1', 'c2', 'c3'],
                              guide_id=['A-1', 'A-1', 'NT-1']))
    labels = pd.DataFrame(dict(rna_class=['KO', 'NP', 'NT']), index=['c1', 'c2', 'c3'])
    result, _ = summarize_guides(catalogue, cells, labels)
    assert result.iloc[0].model_RNA_response_fraction == .5
    assert result.iloc[0].n_model_evaluated == 2
    assert result.iloc[1].status == 'no_QC_records'
    assert result.iloc[3].status == 'NTC_reference_not_scored'
    assert np.isnan(result.iloc[1].model_RNA_response_fraction)
    assert np.isnan(result.iloc[3].model_RNA_response_fraction)


def test_classification_identity_is_checked():
    catalogue = pd.DataFrame(dict(guide=['A-1'], target_gene_name=['A'], source_control_type=['targeting']))
    cells = pd.DataFrame(dict(cell_barcode=['c1'], guide_id=['A-1']))
    with pytest.raises(ValueError, match='axes'):
        summarize_guides(catalogue, cells, pd.DataFrame(dict(rna_class=['KO']), index=['other']))


def test_checkpoint_strings_round_trip(tmp_path):
    import anndata as ad
    a = ad.AnnData(np.eye(3), obs=pd.DataFrame({'guide': pd.array(['A-1', 'A-2', 'NT-1'], dtype='string')}))
    expected = a.obs.guide.tolist()
    checkpoint_strings(a).write_h5ad(tmp_path/'counts.h5ad')
    b = ad.read_h5ad(tmp_path/'counts.h5ad')
    assert b.obs.guide.tolist() == expected
    assert np.array_equal(a.X, b.X)


def test_official_mixscape_interface_with_known_synthetic_response():
    """Integration smoke of the actual published package, not efficacy validation."""
    pt = pytest.importorskip('pertpy')
    import anndata as ad
    import scanpy as sc
    from scipy import sparse

    rng = np.random.default_rng(9128)
    x = rng.poisson(2, size=(180, 256)).astype(np.float32)
    x[90:150, :20] += rng.poisson(20, size=(60, 20))
    a = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame(
        {'group': ['NT'] * 90 + ['synthetic_target'] * 90},
        index=[f'c{i}' for i in range(180)]))
    sc.pp.normalize_total(a, target_sum=10000)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=40, random_state=20261029)
    ms = pt.tl.Mixscape()
    ms.perturbation_signature(a, 'group', 'NT', n_neighbors=20, n_dims=40,
                             random_state=20261029, n_jobs=2)
    ms.mixscape(a, 'group', 'NT', layer='X_pert', random_state=20261029)
    assert a.obs.iloc[:90].mixscape_class_global.eq('NT').all()
    assert a.obs.mixscape_class_global.isin(['NT', 'NP', 'KO']).all()
    assert np.isfinite(a.obs.mixscape_class_p_ko).all()
    assert a.obs.iloc[90:150].mixscape_class_global.eq('KO').any()


def test_interrupted_classification_reuses_saved_input(tmp_path, monkeypatch):
    """An actual official run resumes with source access explicitly forbidden."""
    pytest.importorskip('pertpy')
    from scipy import sparse
    rng = np.random.default_rng(121)
    x = rng.poisson(2, (180, 256)).astype(np.int32)
    x[120:160, :20] += 20
    cells = pd.DataFrame(dict(
        cell_barcode=[f'c{i}' for i in range(180)],
        guide_id=[f'NT-{i//12}' for i in range(120)] + ['A-1']*60,
        total_counts=x.sum(1), role=['non_targeting_control']*120 + ['targeting']*60,
        target=['']*120 + ['A']*60))
    genes = pd.DataFrame(dict(gene_id=[f'g{i}' for i in range(256)],
                              gene_name=[f'g{i}' for i in range(256)]))
    catalogue = pd.DataFrame(dict(guide=[f'NT-{i}' for i in range(10)]+['A-1','absent'],
                                  target_gene_name=['']*10+['A','B'],
                                  source_control_type=['non_targeting_control']*10+['targeting']*2))
    path = tmp_path/'catalogue.csv'
    catalogue.to_csv(path,index=False)
    spec = dict(id='synthetic', input=dict(catalogue=str(path), input_scope='synthetic_only'))
    monkeypatch.setattr(workflow, 'load_counts', lambda *_: (sparse.csr_matrix(x),cells,genes,{'synthetic':True}))
    official = workflow.apply_official
    def stop_after_checkpoint(_):
        raise RuntimeError('intentional classification interruption')
    monkeypatch.setattr(workflow, 'apply_official', stop_after_checkpoint)
    with pytest.raises(RuntimeError, match='intentional'):
        workflow.run_dataset(spec,tmp_path/'result',lambda *_a,**_k:None)
    assert (tmp_path/'result/checkpoint_receipt.json').exists()
    def no_source(*_):
        raise AssertionError('The count source was reread')
    monkeypatch.setattr(workflow, 'load_counts', no_source)
    monkeypatch.setattr(workflow, 'apply_official', official)
    result = workflow.run_dataset(spec,tmp_path/'result',lambda *_a,**_k:None)
    assert result['design_guides']==12
    assert result['DNA_editing_efficiency_assessed'] is False
    assert result['held_whole_NTC_guide_diagnostic']['queried_guides']==2
    assert workflow.run_dataset(spec,tmp_path/'result',lambda *_a,**_k:None)==result
    table=pd.read_csv(tmp_path/'result/all_design_guides_RNA_model.csv')
    assert table.loc[table.guide.eq('absent'),'status'].item()=='no_QC_records'
