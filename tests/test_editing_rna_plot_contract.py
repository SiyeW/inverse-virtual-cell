import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from plot_editing_rna_cells import join_qc_counts, evidence_arrays


def test_each_cohort_uses_its_actual_assigned_record_denominator():
    data={'company':{'guides':[dict(guide='g',assigned_QC_records=31,observed_records=3)]},
          'chip01':{'guides':[dict(guide='g',assigned_QC_records=31,observed_records=2)]}}
    qc=pd.DataFrame([dict(dataset='company',guide='g',records=31),dict(dataset='chip01',guide='g',records=5)])
    changed=join_qc_counts(data,qc)
    assert data['chip01']['guides'][0]['assigned_QC_records']==5 and len(changed)==1


def test_missing_site_observation_is_unknown_not_reference():
    labels,indel,ref=evidence_arrays([],['g'],['cell_without_site_RNA'])
    assert labels[0]=='No site information' and indel[0]==0 and ref[0]==0
