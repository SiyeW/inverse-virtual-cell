from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from plot_guided_rna_cohort import cohort_summary


def example():
    design=[dict(guide=f'G{i}',source_control_type='targeting') for i in range(6)]
    regions=[dict(guides=[dict(guide=f'G{i}') for i in range(4)])]
    c=dict(dataset_id='toy',label='fixture',candidate=dict(all_target_guides=6,observed_target_guides=5,uniquely_placed_observed=4,regions=regions),
           guide_rows=[dict(guide=f'G{i}',own_RNA_families_observed=0 if i in [0,3] else 3,scan_capped=i==0) for i in range(4)],variant_rows=[])
    return design,c


def test_missing_capped_and_unmapped_are_not_negative_edit_calls():
    d,c=example();s=cohort_summary(d,c)
    assert s['counts']==[1,1,1,1,2,0,0]
    assert sum(s['counts'])==6 and not s['functional_KO_efficiency_assessed']


def test_deletion_support_counts_even_when_whole_core_deleted():
    d,c=example()
    c['variant_rows']=[dict(guide='G3',kind='DEL',guide_alt=2,guide_ref=0,guide_other=0,NTC_alt=0,NTC_ref=0,NTC_other=0)]
    s=cohort_summary(d,c)
    assert s['counts'][6]==1 and s['counts'][3]==0
    assert s['indel_rows_without_NTC_opportunity']==1 and s['points']==[]


def test_incomplete_cohort_cannot_be_drawn_as_whole():
    d,c=example();c['guide_rows'].pop()
    with pytest.raises(ValueError,match='Incomplete'):cohort_summary(d,c)
