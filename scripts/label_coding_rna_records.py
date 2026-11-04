"""Join project-ORF consequences to existing per-barcode RNA allele evidence.

All available guides under the original 104 targets are retained, separately in
company merged records and the independent first capture. No DNA genotype or
functional-KO label is assigned. Stop/start-region changes are kept separately
from internal frameshift evidence. No BAM is opened.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from annotate_rna_indel_orfs import classify, references
from explore_editing_rna_cells import add_same_indel_support

MODES = ("internal_frameshift", "two_distinct_UMIs", "not_observed_in_NTC")


def key(row):
    return (row["guide"], row["chrom"], int(row["pos0"]), row["ref"], row["alt"])


def annotate(cohort: dict, models: dict):
    annotated = []
    identities = set()
    for variant in cohort["variants"]:
        if key(variant) in identities:
            raise ValueError("Duplicate accepted guide-variant identity")
        identities.add(key(variant))
        row = {k: variant[k] for k in (
            "guide", "target_gene_id", "chrom", "pos0", "ref", "alt", "kind", "length_change",
            "NTC_alt", "guide_alt", "cell_support", "primary_197")}
        if row["target_gene_id"] not in models:
            raise ValueError("Variant target lacks the accepted project ORF reference")
        if row["kind"] in ("INS", "DEL"):
            row.update(classify(row, models[row["target_gene_id"]]))
        else:
            row.update(category="Substitution", reason="not_an_indel", start_or_stop_overlap=False)
        row["internal_frameshift"] = row["category"] == "Frameshift" and not row["start_or_stop_overlap"]
        annotated.append(row)
    return annotated


def label_records(cohort: dict, annotated: list[dict]):
    cells = [dict(r) for r in cohort["cells"]]
    if len({(r["guide"],r["cell_barcode"]) for r in cells}) != len(cells):
        raise ValueError("Repeated guide/barcode evidence row")
    internal = [r for r in annotated if r["internal_frameshift"]]
    # Reuse the accepted conservative bound for one CB/UB with multiple GXs.
    focused = add_same_indel_support(dict(cells=cells, guides=[], display=cohort["display"], variants=internal))
    ref_max, ref_clean_max, alt_clean_max = defaultdict(int), defaultdict(int), defaultdict(int)
    for variant in internal:
        for cell in variant["cell_support"]:
            k = (variant["guide"],cell["cell_barcode"])
            ref_max[k] = max(ref_max[k],cell["ref"])
            if not variant["NTC_alt"]:
                ref_clean_max[k] = max(ref_clean_max[k],cell["ref"])
                alt_clean_max[k] = max(alt_clean_max[k],cell["alt"])
    # Estimate the same conservative duplication allowance for reference families.
    duplicates = defaultdict(int)
    for guide in cohort["display"]:
        molecular_genes = defaultdict(set)
        for molecule in guide["molecules"]:
            molecular_genes[(molecule["cell_barcode"],molecule["UMI"])].add(molecule["GX"])
        for (barcode,_),genes in molecular_genes.items():
            duplicates[(guide["guide"],barcode)] += max(0,len(genes)-1)
    output = []
    for cell in focused["cells"]:
        k = (cell["guide"],cell["cell_barcode"])
        row = {name:cell[name] for name in (
            "guide", "cell_barcode", "target_gene_id", "regional_families", "cut_ref_families",
            "indel_families", "SNV_families", "ambiguous_families")}
        row.update(frameshift_max_same_variant_families=cell["max_same_indel_families"],
                   frameshift_distinct_UMI_lower_bound=cell["same_indel_distinct_UMI_lower_bound"],
                   frameshift_site_reference_max=ref_max[k],
                   frameshift_site_reference_UMI_lower_bound=max(0,ref_max[k]-duplicates[k]),
                   NTC_unobserved_frameshift_max=alt_clean_max[k],
                   NTC_unobserved_site_reference_max=ref_clean_max[k])
        row["display_category"] = display_category(row)
        output.append(row)
    return output


def display_category(row):
    if row.get("frameshift_max_same_variant_families",0) > 0:
        return "Internal frameshift RNA"
    if row.get("indel_families",0) > 0 or row.get("SNV_families",0) > 0:
        return "Other variant RNA"
    if row.get("frameshift_site_reference_max",0) > 0 and not row.get("ambiguous_families",0):
        return "Reference RNA at frameshift sites"
    return "Insufficient site information"


def membership(row, mode):
    if mode not in MODES:
        raise ValueError("Unknown comparison mode")
    # Conflicting families remain in the figures; they do not enter comparisons.
    if row.get("ambiguous_families",0):
        return False,False
    if mode == "two_distinct_UMIs":
        positive = row.get("frameshift_distinct_UMI_lower_bound",0) >= 2
        reference = row.get("frameshift_site_reference_UMI_lower_bound",0) >= 2
    elif mode == "not_observed_in_NTC":
        positive = row.get("NTC_unobserved_frameshift_max",0) > 0
        reference = row.get("NTC_unobserved_site_reference_max",0) > 0
    else:
        positive = row.get("frameshift_max_same_variant_families",0) > 0
        reference = row.get("frameshift_site_reference_max",0) > 0
    negative = bool(reference and not row.get("indel_families",0) and not row.get("SNV_families",0))
    if positive and negative:
        raise ValueError("A record cannot enter both evidence classes")
    return bool(positive),negative


def build(sequence, reference, previous):
    models = references(reference)
    if len(models) != 104:
        raise ValueError("The complete accepted 104-target reference is required")
    result = {"scope":"All available guides and barcode records under the original 104 targets; cohort-specific RNA evidence, not DNA genotypes or functional KO calls.",
              "target_ids":sorted(models), "modes":MODES,
              "reference_definition":"Reference RNA at observed internal-frameshift candidate sites for the SAME guide; no observed indel/substitution/conflict in that record. It is not proven DNA wild type.",
              "two_UMI_sensitivity":"Both positive and reference evidence require the accepted conservative >=2 distinct-UMI lower bound at the same candidate site.",
              "start_stop_changes":"Retained as other variant RNA; not used as internal frameshift evidence."}
    for dataset in ("company","chip01"):
        annotated = annotate(sequence[dataset],models)
        records = label_records(sequence[dataset],annotated)
        result[dataset] = {"variants":annotated,"records":records,
                           "variant_categories":dict(Counter(r["category"] for r in annotated)),
                           "internal_frameshift_variant_observations":sum(r["internal_frameshift"] for r in annotated)}
    old = {key(r):r for r in previous["events"]}
    primary = {key(r):r for r in result["company"]["variants"] if r["primary_197"] and r["kind"] in ("INS","DEL")}
    if set(primary) != set(old) or any(primary[k]["category"] != old[k]["category"] for k in old):
        raise ValueError("Primary 145-event coding predictions changed during the record join")
    result["primary_145_event_predictions_unchanged"] = True
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("sequences","references","previous","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    parser.add_argument("--commit",required=True)
    args=parser.parse_args();load=lambda p:json.loads(p.read_text(encoding="utf-8-sig"))
    result=build(load(args.sequences),load(args.references),load(args.previous))
    result["code_commit"]=args.commit
    result["input_sha256"]={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.sequences,args.references,args.previous)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists():raise FileExistsError("Preserve previous label packets")
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({d:{k:v for k,v in result[d].items() if k not in ("variants","records")} for d in ("company","chip01")}))


if __name__=="__main__":main()
