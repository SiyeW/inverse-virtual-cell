import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_cluster_reference_coverage import analyze


def fixture():
    cells=pd.DataFrame({'guide_id':['n1','n1','n2','n2','a1','a1']},index=list('abcdef'))
    records=pd.DataFrame(dict(guide=['n1','n2','a1'],cells=[2,2,2],source_control_type=['non_targeting_control']*2+['targeting']))
    clusters=pd.DataFrame(dict(Barcode=list('fedcba'),Cluster=['2','1','2','1','2','1']))
    return cells,records,clusters


def test_label_join_uses_barcodes_and_conserves_counts():
    cells,records,clusters=fixture();table,guides,result=analyze(cells,records,clusters)
    assert result['matched_analysis_cells']==6 and result['ntc_cells']==4
    assert result['clusters_with_ntc']==2
    assert table.total_cells.sum()==6
    np.testing.assert_allclose(guides.expression_cluster_total_variation,0)
    assert guides.reference_overlap_cells.eq(0).all()


def test_missing_or_duplicate_labels_are_not_silently_dropped():
    cells,records,clusters=fixture()
    with pytest.raises(ValueError,match='Missing'):analyze(cells,records,clusters.iloc[:-1])
    with pytest.raises(ValueError,match='Nonunique'):analyze(cells,records,pd.concat([clusters,clusters.iloc[:1]]))


def test_report_exposes_a_cluster_without_ntc_without_inventing_a_cell_type():
    cells,records,clusters=fixture();clusters.loc[clusters.Barcode.isin(['e','f']),'Cluster']='3'
    table,guides,result=analyze(cells,records,clusters)
    assert result['clusters_with_ntc']==2 and result['expression_clusters']==3
    assert result['cluster_ntc_fraction_range'][0]==0
    assert not result['clusters_are_validated_cell_types'] and not result['clusters_are_capture_labels']


def test_old_catalogue_counts_are_not_reused():
    cells,records,clusters=fixture();records['cells']=999999
    records['selected_as_stable_control']=True
    table,guides,result=analyze(cells,records,clusters)
    assert result['analysis_cells']==6 and result['ntc_cells']==4
    assert guides.cells.sum()==6 and result['guide_counts_recomputed_from_analysis_cells']
