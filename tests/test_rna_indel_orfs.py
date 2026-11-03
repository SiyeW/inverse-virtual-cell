import importlib.util
from pathlib import Path
import sys


path = Path(__file__).parents[1] / "scripts" / "annotate_rna_indel_orfs.py"
spec = importlib.util.spec_from_file_location("rna_orf_consequences", path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def event(pos, ref, alt, chrom="chrF"):
    return dict(pos0=pos, ref=ref, alt=alt, chrom=chrom, length_change=len(alt) - len(ref))


def positive():
    return module.ORF("chrF", "+", ((100, 115),), "ATGAAACCCGGGTAA")


def test_vcf_anchor_outside_cds_does_not_hide_deleted_coding_base():
    call = module.classify(event(99, "CAT", "C"), [positive()])
    assert call["category"] == "Frameshift"
    assert call["start_or_stop_overlap"]


def test_in_frame_insertion_and_frameshift_insertion():
    assert module.classify(event(105, "A", "AAAA"), [positive()])["category"] == "In-frame"
    assert module.classify(event(105, "A", "AA"), [positive()])["category"] == "Frameshift"


def test_right_anchored_insertion_uses_its_actual_reference_position():
    assert module.classify(event(106, "C", "TC"), [positive()])["category"] == "Frameshift"


def test_negative_strand_reference_is_compared_in_genomic_orientation():
    model = module.ORF("chrF", "-", ((100, 115),), "ATGAAACCCGGGTAA")
    assert module.classify(event(103, "CCC", "C"), [model])["category"] == "Frameshift"
    assert module.classify(event(102, "ACCC", "A"), [model])["category"] == "In-frame"


def test_deletion_partly_outside_cds_is_not_classified_from_total_length():
    call = module.classify(event(97, "CCCAT", "C"), [positive()])
    assert call["category"] == "Unresolved"
    assert "boundary" in call["reason"]


def test_outer_boundary_insertion_is_not_a_coding_frameshift():
    assert module.classify(event(99, "C", "CA"), [positive()])["category"] == "Outside ORF"
    assert module.classify(event(114, "A", "AT"), [positive()])["category"] == "Outside ORF"


def test_internal_splice_boundary_is_not_treated_as_contiguous_coding():
    model = module.ORF("chrF", "+", ((100, 106), (200, 209)), "ATGAAACCCGGGTAA")
    call = module.classify(event(105, "A", "AA"), [model])
    assert call["category"] == "Unresolved"
    assert "splice" in call["reason"]


def test_multiple_references_require_agreement_instead_of_best_result():
    different = module.ORF("chrF", "+", ((200, 215),), "ATGAAACCCGGGTAA")
    assert module.classify(event(105, "A", "AA"), [positive(), different])["category"] == "Unresolved"
    assert module.classify(event(105, "A", "AA"), [positive(), positive()])["category"] == "Frameshift"


def test_missing_reference_and_ref_mismatch_stay_unresolved():
    assert module.classify(event(105, "A", "AA"), [])["category"] == "Unresolved"
    assert module.classify(event(105, "G", "GG"), [positive()])["category"] == "Unresolved"
