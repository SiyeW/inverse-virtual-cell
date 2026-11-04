"""Build sequence-checked references from the project's intended ORF annotation.

The preparation step writes a small indexed-FASTA request; validation consumes
only that returned packet. It never opens a BAM or an RNA expression matrix.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

from plot_design_cut_positions import (
    Model, boundary_offset, merge_blocks, read_models, reverse_complement,
    supplementary_model,
)

STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "CTG", "GTG", "TTG", "ACG", "ATT", "ATC", "ATA", "AGG", "AAG"}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def slice_transcript(exons, strand, start, end):
    """Return genomic blocks for a 0-based half-open spliced-RNA interval."""
    if strand not in ("+", "-") or not 0 <= start < end:
        raise ValueError("Invalid transcript interval")
    prefix, selected = 0, []
    for left, right in exons if strand == "+" else reversed(exons):
        length = right-left
        low, high = max(start-prefix, 0), min(end-prefix, length)
        if low < high:
            selected.append((left+low, left+high) if strand == "+" else (right-high, right-low))
        prefix += length
    if sum(e-s for s,e in selected) != end-start:
        raise ValueError("Transcript interval exceeds annotated exons")
    return merge_blocks(selected)


def following_codon(model):
    boundary = model.coding[-1][1] if model.strand == "+" else model.coding[0][0]
    offset = boundary_offset(model.exons, model.strand, boundary)
    if offset is None:
        raise ValueError("Coding end is not an exonic boundary")
    return slice_transcript(model.exons, model.strand, offset, offset+3)


def collect_request(design, primary, gene_table, supplements, parents):
    cohort = primary["records"]
    if len(cohort) != len({r["guide"] for r in cohort}) or len(cohort) != 197:
        raise ValueError("Expected the complete accepted primary 197 guides")
    ids = {r["guide"] for r in cohort}
    chosen = [r for r in design["records"] if r["guide"] in ids]
    if len(chosen) != 197 or any(r["status"] != "assessed" for r in chosen):
        raise ValueError("Incomplete coordinate validation for primary cohort")
    by_name = dict(gene_table)
    supplemental = {}
    for row in supplements["models"]:
        _, model = supplementary_model(row, parents["current_transcripts"])
        if model.name in by_name:
            raise ValueError("Supplement would replace a project model")
        by_name[model.name] = model
        supplemental[model.name] = row
    per_gene = defaultdict(set)
    for row in chosen:
        per_gene[row["target_gene_id"]].update(row["annotation_models"])
    models, regions = [], set()
    for gene, names in sorted(per_gene.items()):
        # The supplied project annotation defines one intended model per target.
        # An accidental ambiguity must be investigated, not silently selected.
        if len(names) != 1:
            raise ValueError(f"Intended project ORF identity is ambiguous: {gene}")
        model = by_name[next(iter(names))]
        blocks = model.coding
        expected_start = model.name.rsplit("|", 1)[1] if model.name.startswith("PB.") else None
        source_kind, omitted_stop = "project_genePred", False
        if model.name in supplemental:
            raw = supplemental[model.name]
            source_kind = raw["format"]
            expected_start = "ATG"
            if source_kind == "Ensembl_expanded_transcript":
                omitted_stop = model.coding_length == raw["raw"]["Translation"]["length"]*3
        extension = following_codon(model) if omitted_stop else ()
        row = {"target_gene_id":gene,"source_model_name":model.name,"source_kind":source_kind,
               "chrom":model.chrom,"strand":model.strand,"blocks":list(blocks),
               "expected_start_codon":expected_start,"source_coding_length":model.coding_length,
               "source_omits_terminal_stop":omitted_stop,"terminal_stop_blocks":list(extension)}
        models.append(row)
        regions.update((model.chrom,s,e) for s,e in (*blocks,*extension))
    targets = []
    for row in cohort:
        strand = row["original_spacer_strand"]
        left, right = row["start0"], row["end0"]
        interval = (left,right+3) if strand == "+" else (left-3,right)
        if strand not in ("+", "-") or interval[1]-interval[0] != 23:
            raise ValueError("Invalid original 23nt guide target")
        targets.append({"guide":row["guide"],"chrom":row["chrom"],"start0":interval[0],
                        "end0":interval[1],"strand":strand,"sequence":row["original_Target_sequence"]})
        regions.add((row["chrom"],*interval))
    return {"models":models,"original_targets":targets,"zero_based_regions":sorted(regions),
            "reference_bases_requested":sum(e-s for _,s,e in regions),
            "scope":"All 104 intended targets of the original 197 sgRNAs; original project annotation plus exact positive-control transcripts. No response-based model choice."}


def read_fasta(path):
    sequences, name, pieces = {}, None, []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if name is not None:
                if name in sequences:raise ValueError("Repeated reference region")
                sequences[name] = "".join(pieces).upper()
            name, pieces = line[1:].strip(), []
        else:
            pieces.append(line.strip())
    if name is not None:
        if name in sequences:raise ValueError("Repeated reference region")
        sequences[name] = "".join(pieces).upper()
    return sequences


def label(chrom,start,end):
    return f"{chrom}:{start+1}-{end}"


def oriented_sequence(chrom, strand, blocks, sequences):
    sequence = "".join(sequences[label(chrom,s,e)] for s,e in blocks)
    return sequence if strand == "+" else reverse_complement(sequence)


def validate_reference(row, sequences):
    blocks = tuple(map(tuple,row["blocks"]))
    sequence = oriented_sequence(row["chrom"],row["strand"],blocks,sequences)
    if row["source_omits_terminal_stop"]:
        extension = tuple(map(tuple,row["terminal_stop_blocks"]))
        sequence += oriented_sequence(row["chrom"],row["strand"],extension,sequences)
        blocks = merge_blocks((*blocks,*extension))
    errors = []
    if len(sequence)%3:
        errors.append("coding_length_not_multiple_of_three")
    if set(sequence)-set("ACGT"):
        errors.append("ambiguous_reference_bases")
    expected = row["expected_start_codon"]
    if expected and sequence[:3] != expected:
        errors.append("reference_start_differs_from_original_annotation")
    if not expected and sequence[:3] not in STARTS:
        errors.append("unexpected_reference_initiation_codon")
    if sequence[-3:] not in STOPS:
        errors.append("missing_annotated_terminal_stop")
    internal = [i for i in range(3,len(sequence)-3,3) if sequence[i:i+3] in STOPS]
    if internal:
        errors.append("internal_in_frame_reference_stop")
    return dict(row,validated_blocks=blocks,reference_ORF_sequence=sequence,
                reference_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                reference_ORF_length=len(sequence),reference_start_codon=sequence[:3],
                reference_terminal_codon=sequence[-3:],internal_stop_offsets=internal,
                sequence_check_errors=errors,reference_sequence_consistent_with_design=not errors)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="stage",required=True)
    prepare=sub.add_parser("request")
    for name in ["design_positions","primary_join","gene_table","supplements","parents"]:
        prepare.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    prepare.add_argument("--output",type=Path,required=True)
    prepare.add_argument("--commit",required=True)
    accept=sub.add_parser("validate")
    for name in ["request","fasta","extraction_receipt","output"]:
        accept.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    accept.add_argument("--commit",required=True)
    args=parser.parse_args()
    load=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
    if args.stage=="request":
        design,primary=load(args.design_positions),load(args.primary_join)
        if not design["original_relative_positions_equal_for_all_assessed"]:
            raise ValueError("Unaccepted original cut-site positions")
        for name in ["gene_table","supplementary_models","enst_parents"]:
            path={"gene_table":args.gene_table,"supplementary_models":args.supplements,"enst_parents":args.parents}[name]
            if digest(path)!=design["source_files"][name]["sha256"]:
                raise ValueError("Source annotation changed after design acceptance")
        request=collect_request(design,primary,read_models(args.gene_table),load(args.supplements),load(args.parents))
        if request["reference_bases_requested"]>2*1024*1024:
            raise ValueError("Unexpectedly large indexed reference request")
        request["code_commit"]=args.commit
        request["input_sha256"]={str(p):digest(p) for p in [args.design_positions,args.primary_join,args.gene_table,args.supplements,args.parents]}
        args.output.mkdir(parents=True,exist_ok=False)
        (args.output/"request.json").write_text(json.dumps(request,ensure_ascii=False,indent=2),encoding="utf-8")
        (args.output/"regions.txt").write_text("".join(label(*r)+"\n" for r in request["zero_based_regions"]),encoding="ascii")
        print(json.dumps({"targets":len(request["models"]),"guides":len(request["original_targets"]),"regions":len(request["zero_based_regions"]),"bases":request["reference_bases_requested"]}))
    else:
        request,receipt=load(args.request),load(args.extraction_receipt)
        if digest(args.fasta)!=receipt["output_sha256"]:
            raise ValueError("Extracted FASTA differs from accepted source SHA")
        if receipt["request_sha256"]!=digest(args.request) or not receipt["source_identity_unchanged"]:
            raise ValueError("Unaccepted indexed extraction identity")
        sequences=read_fasta(args.fasta)
        if set(sequences)!={label(*r) for r in request["zero_based_regions"]}:
            raise ValueError("Missing or extra extracted reference regions")
        for chrom,s,e in request["zero_based_regions"]:
            if len(sequences[label(chrom,s,e)])!=e-s:
                raise ValueError("Extracted interval length mismatch")
        for target in request["original_targets"]:
            actual=sequences[label(target["chrom"],target["start0"],target["end0"])]
            if target["strand"]=="-":actual=reverse_complement(actual)
            if actual!=target["sequence"]:
                raise ValueError("Original guide target no longer matches the reference")
        models=[validate_reference(row,sequences) for row in request["models"]]
        accepted=[row for row in models if row["reference_sequence_consistent_with_design"]]
        result={"code_commit":args.commit,"all_197_original_23nt_targets_match_GRCh38":True,
                "tested_models":len(models),"sequence_consistent_models":len(accepted),
                "accepted_target_count":len({r["target_gene_id"] for r in accepted}),
                "scope":request["scope"],"models":models,
                "source_reference_receipt":str(args.extraction_receipt),
                "input_sha256":{str(p):digest(p) for p in [args.request,args.fasta,args.extraction_receipt]}}
        args.output.parent.mkdir(parents=True,exist_ok=True)
        if args.output.exists():raise FileExistsError("Preserve accepted reference packets")
        args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
        print(json.dumps({k:v for k,v in result.items() if k not in ("models","input_sha256")},ensure_ascii=False))


if __name__=="__main__":main()
