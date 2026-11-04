"""Coordinate and identity tests for intended coding cut-site distances."""
import importlib.util
from pathlib import Path
import sys

import pytest

spec = importlib.util.spec_from_file_location("cut_positions",Path(__file__).parents[1]/"scripts/plot_design_cut_positions.py")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_plus_excludes_intron_and_counts_boundary_not_base_index():
    assert module.boundary_offset(((100,110),(200,220)),"+",201) == 11


def test_minus_orders_exons_in_transcript_direction():
    assert module.boundary_offset(((100,110),(200,220)),"-",109) == 21


def test_intronic_cut_not_mapped_to_nearest_coding_exon():
    assert module.boundary_offset(((100,110),(200,220)),"+",150) is None


@pytest.mark.parametrize("strand,start,stop",[("+",100,220),("-",220,100)])
def test_outer_boundaries_retained(strand,start,stop):
    blocks=((100,110),(200,220))
    assert module.boundary_offset(blocks,strand,start) == 0
    assert module.boundary_offset(blocks,strand,stop) == 30


def test_gtf_start_is_converted_once_and_stop_is_added(tmp_path):
    model=module.Model("A","chr1","+",((100,110),),((101,110),))
    gtf=tmp_path/"source.gtf"
    attributes='gene_id "A"; transcript_id "A";'
    gtf.write_text("\n".join(f'chr1\tsrc\t{kind}\t{start}\t{end}\t.\t+\t.\t{attributes}'
                            for kind,start,end in [("exon",101,110),("CDS",102,107),
                                                   ("start_codon",102,104),("stop_codon",108,110)])+"\n")
    assert module.check_gtf(gtf,{"A":model})["models_checked"]==1
    gtf.write_text(gtf.read_text().replace("\texon\t101\t","\texon\t100\t"))
    with pytest.raises(ValueError,match="exon mismatch"):
        module.check_gtf(gtf,{"A":model})


def test_pb_name_preserves_transcript_identity_but_not_feature_omitted_strand():
    assert module.identity("PB.27.1:chr1:+|9|916:110:311|uORF|ATG") == "PB.27.1_chr1_9_916_110_311_uORF_ATG"


def test_duplicate_alias_is_not_a_second_guide_and_source_distance_is_checked():
    feature={"guide":"TYMSOS-1","sequence":"ACGT","target_gene_id":"ENSG1",
             "source_control_type":"targeting","single_guide_after_stage1_qc":"0"}
    base={"chrom":"chr1","cut0":103,"relative":3,"original_row":1}
    rows={"ACGT":[dict(base,original_gene="ENSG1"),dict(base,original_gene="PB.alias",original_row=2)]}
    models={"ENSG1":[module.Model("ENSG1","chr1","+",((100,112),),((100,112),))]}
    actual=module.join_feature(feature,rows,models,{})
    assert actual["distance"] == 3 and actual["RNA_QC_records"] == 0
    assert actual["source_rows_with_this_spacer"] == [1,2]
    rows["ACGT"][0]["relative"]=4
    with pytest.raises(ValueError,match="disagrees"):
        module.join_feature(feature,rows,models,{})


def test_source_transcript_parent_is_an_identity_join_only():
    feature={"guide":"PC-1","sequence":"ACGT","target_gene_id":"ENSG1",
             "source_control_type":"positive_control","single_guide_after_stage1_qc":"12"}
    rows={"ACGT":[{"original_gene":"ENST1.2","chrom":"chr1","cut0":103,"relative":3,"original_row":1}]}
    models={"ENSG1":[module.Model("ENSG1","chr1","+",((100,112),),((100,112),))]}
    assert module.join_feature(feature,rows,models,{"ENST1":{"Parent":"ENSG1"}})["distance"]==3
