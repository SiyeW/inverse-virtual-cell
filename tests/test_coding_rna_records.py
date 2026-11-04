from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).parents[1]/"scripts"))
from label_coding_rna_records import annotate,label_records,membership,display_category
from annotate_rna_indel_orfs import ORF


def record(**kwargs):
    row=dict(guide="g",cell_barcode="c",target_gene_id="t",regional_families=3,
             cut_ref_families=1,indel_families=0,SNV_families=0,ambiguous_families=0)
    row.update(kwargs);return row


def variant(pos,alt,cells,ntc=0):
    return dict(guide="g",target_gene_id="t",chrom="chrF",pos0=pos,ref="A",alt=alt,
                kind="INS",length_change=len(alt)-1,NTC_alt=ntc,guide_alt=sum(r["alt"] for r in cells),
                cell_support=cells,primary_197=True)


def support(cb,alt=0,ref=0):
    return dict(cell_barcode=cb,alt=alt,ref=ref,other=0,ambiguous=0)


def test_termination_frameshift_retained_but_not_an_internal_frameshift_label():
    models={"t":[ORF("chrF","+",((100,112),),"ATGAAACCCTAA")]}
    a=annotate(dict(variants=[variant(105,"AA",[support("c",alt=1)]),
                             variant(110,"AA",[support("d",alt=1)])]),models)
    assert a[0]["internal_frameshift"] and a[0]["category"]=="Frameshift"
    assert not a[1]["internal_frameshift"] and a[1]["start_or_stop_overlap"]


def test_reference_requires_coverage_at_the_matching_guide_variant_site():
    v=dict(kind="INS",internal_frameshift=True,guide="g",NTC_alt=0,cell_support=[support("p",alt=1),support("r",ref=1)])
    cohort=dict(cells=[record(cell_barcode="p",indel_families=1),record(cell_barcode="r"),record(cell_barcode="window_only")],
                display=[],variants=[])
    labels={r["cell_barcode"]:r for r in label_records(cohort,[v])}
    assert membership(labels["p"],"internal_frameshift")== (True,False)
    assert membership(labels["r"],"internal_frameshift")== (False,True)
    assert membership(labels["window_only"],"internal_frameshift")== (False,False)
    assert display_category(labels["window_only"])=="Insufficient site information"


def test_two_UMI_sensitivity_applies_to_reference_and_variant_evidence():
    x=record(frameshift_distinct_UMI_lower_bound=0,frameshift_site_reference_UMI_lower_bound=1)
    assert membership(x,"two_distinct_UMIs")== (False,False)
    x["frameshift_site_reference_UMI_lower_bound"]=2
    assert membership(x,"two_distinct_UMIs")== (False,True)


def test_conflicting_and_other_variant_records_do_not_become_reference_controls():
    x=record(frameshift_site_reference_max=2,ambiguous_families=1)
    assert membership(x,"internal_frameshift")== (False,False)
    x.update(ambiguous_families=0,SNV_families=1)
    assert membership(x,"internal_frameshift")== (False,False)
    assert display_category(x)=="Other variant RNA"


def test_ntc_sensitivity_does_not_turn_shared_variant_support_into_a_negative():
    x=record(indel_families=1,frameshift_max_same_variant_families=2,
             NTC_unobserved_frameshift_max=0,NTC_unobserved_site_reference_max=1)
    assert membership(x,"internal_frameshift")== (True,False)
    assert membership(x,"not_observed_in_NTC")== (False,False)


def test_same_UMI_in_two_gene_tags_is_not_two_distinct_RNA_UMIs():
    v=dict(kind="INS",internal_frameshift=True,guide="g",NTC_alt=0,cell_support=[support("c",alt=2)])
    cohort=dict(cells=[record(indel_families=2)],display=[dict(guide="g",molecules=[
        dict(cell_barcode="c",UMI="U",GX="A"),dict(cell_barcode="c",UMI="U",GX="B")])])
    x=label_records(cohort,[v])[0]
    assert x["frameshift_max_same_variant_families"]==2
    assert x["frameshift_distinct_UMI_lower_bound"]==1
    assert membership(x,"two_distinct_UMIs")== (False,False)
