import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from joint_target_rna_region import combine,check_membership,summarize

def fixture():
    e=pd.DataFrame(dict(guide=list('abcdef'),target_gene_id=list('ABCDEF'),target=list('ABCDEF'),cells=[0,1,1,1,1,1],
        detected_cells=[0,1,0,1,0,1],mapping_status=['exact_catalogue_id']*6,reference_status=['not_detected']*6))
    r=pd.DataFrame(dict(guide=list('bcde'),own_RNA_families_observed=[0,2,0,0],max_alt_families=[0,0,0,1],scan_capped=[True,True,False,False]))
    m=pd.DataFrame(dict(guide=list('abcdef'),mapping_status=['unique_exact_genomic_placement']*5+['ambiguous_exact_placement']))
    c=pd.DataFrame(dict(cell_barcode=list('12345'),guide_id=list('bcdef'),total_counts=[20]*5))
    return e,r,m,c

def test_zero_assignment_is_not_RNA_absence_caps_are_not_negative_and_layers_need_not_nest():
    e,r,m,c=fixture();t=combine(e,r,m,c).set_index('guide')
    assert t.loc['a','regional_RNA_status']=='no_assigned_records' and pd.isna(t.loc['a','target_RNA_detected'])
    assert t.loc['b','target_RNA_detected'] and t.loc['b','regional_RNA_status']=='not_assessable' # NTC absent does not exclude it
    assert not t.loc['c','target_RNA_detected'] and t.loc['c','regional_RNA_status']=='observed' and t.loc['c','positive_lower_bound_only']
    assert t.loc['d','regional_RNA_status']=='not_observed'
    assert t.loc['e','regional_RNA_status']=='observed' # complete deletion support is not zero evidence
    assert t.loc['f','regional_RNA_status']=='not_assessable'
    assert sum(x['guides'] for x in summarize(t))==5

def test_identity_record_counts_and_unmapped_features_fail_closed():
    e,r,m,c=fixture();bad=e.copy();bad.loc[1,'cells']=2
    with pytest.raises(ValueError):combine(bad,r,m,c)
    bad=e.copy();bad.loc[1,'mapping_status']='not_found'
    with pytest.raises(ValueError):combine(bad,r,m,c)
    with pytest.raises(ValueError):combine(e,r.iloc[:-1],m,c)
    assert check_membership(c,c.iloc[::-1])==5
    bad=c.copy();bad.loc[0,'guide_id']='z'
    with pytest.raises(ValueError):check_membership(c,bad)
