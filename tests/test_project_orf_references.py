import importlib.util
from pathlib import Path
import sys

import pytest

SCRIPTS=Path(__file__).parents[1]/"scripts"
sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location("project_orf_refs",SCRIPTS/"prepare_project_orf_references.py")
module=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=module
spec.loader.exec_module(module)


def test_transcript_slice_crosses_introns_in_both_directions():
    exons=((100,104),(200,205))
    assert module.slice_transcript(exons,"+",3,6)==((103,104),(200,202))
    assert module.slice_transcript(exons,"-",3,6)==((103,104),(200,202))
    assert module.slice_transcript(exons,"-",0,3)==((202,205),)


def test_following_stop_can_span_two_exons_and_preserves_strand():
    positive=module.Model("p","chrF","+",((100,106),(200,205)),((100,105),))
    assert module.following_codon(positive)==((105,106),(200,202))
    negative=module.Model("n","chrF","-",((100,105),(200,206)),((201,206),))
    assert module.following_codon(negative)==((103,105),(200,201))


def row(strand="+",omitted=False):
    return {"target_gene_id":"t","source_model_name":"source","source_kind":"project_genePred",
            "chrom":"chrF","strand":strand,"blocks":[[100,109]],
            "expected_start_codon":"ATG","source_coding_length":9,
            "source_omits_terminal_stop":omitted,"terminal_stop_blocks":[[109,112]] if omitted else []}


def test_stop_omission_is_only_extended_when_declared_by_source():
    actual=module.validate_reference(row(omitted=True),{"chrF:101-109":"ATGAAACCC","chrF:110-112":"TAA"})
    assert actual["validated_blocks"]==((100,112),)
    assert actual["reference_ORF_sequence"]=="ATGAAACCCTAA"
    assert actual["reference_sequence_consistent_with_design"]
    actual=module.validate_reference(row(),{"chrF:101-109":"ATGAAACCC"})
    assert "missing_annotated_terminal_stop" in actual["sequence_check_errors"]


def test_negative_strand_reference_is_oriented_before_codon_checks():
    actual=module.validate_reference(row("-"),{"chrF:101-109":"TTATTTCAT"})
    assert actual["reference_ORF_sequence"]=="ATGAAATAA"
    assert actual["reference_sequence_consistent_with_design"]


def test_internal_stop_and_inconsistent_source_start_remain_unaccepted():
    r=row();r["blocks"]=[[100,112]];r["source_coding_length"]=12
    actual=module.validate_reference(r,{"chrF:101-112":"CTGTAACCCTAA"})
    assert not actual["reference_sequence_consistent_with_design"]
    assert "internal_in_frame_reference_stop" in actual["sequence_check_errors"]
    assert "reference_start_differs_from_original_annotation" in actual["sequence_check_errors"]


def test_bad_spliced_end_cannot_fabricate_a_stop_extension():
    with pytest.raises(ValueError,match="exceeds"):
        module.slice_transcript(((100,110),),"+",9,12)


def test_duplicate_fasta_headers_are_rejected(tmp_path):
    p=tmp_path/"ref.fa";p.write_text(">chrF:1-3\nATG\n>chrF:1-3\nATG\n")
    with pytest.raises(ValueError,match="Repeated"):
        module.read_fasta(p)
