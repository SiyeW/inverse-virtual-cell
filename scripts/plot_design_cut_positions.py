"""Audit library cut positions against project genePred/GTF coding models.

One observation per unique targeting feature, without RNA-response filtering.
Distances are boundary offsets along spliced coding sequence, not genomic spans.
This evaluates intended design positions and does not measure observed editing.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import statistics


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def identity(name: str) -> str:
    """Join the original ORF name to its strand-free expression feature ID."""
    name = re.sub(r"[:|]", "_", name)
    return re.sub(r"_([+-])_", "_", name) if name.startswith("PB.") else name


def merge_blocks(blocks):
    merged = []
    for start, end in sorted(set(blocks)):
        if not 0 <= start < end:
            raise ValueError("Invalid half-open interval")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def boundary_offset(blocks, strand: str, cut0: int) -> int | None:
    """Count coding bases preceding a genomic cut boundary in transcript order.

    cut0 is a BETWEEN-base coordinate. E.g. a cut after the first base of a
    plus-strand block [100, 110) is cut0=101, and has coding offset 1.
    Splice boundaries have defined offsets, but an intronic cut is not folded
    onto the nearest exon. Both outer ORF boundaries are retained explicitly.
    """
    if strand not in ("+", "-"):
        raise ValueError("Unknown coding strand")
    if not any(start <= cut0 <= end for start, end in blocks):
        return None
    if strand == "+":
        return sum(max(0, min(cut0, end) - start) for start, end in blocks)
    return sum(max(0, end - max(cut0, start)) for start, end in blocks)


@dataclass(frozen=True)
class Model:
    name: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    coding: tuple[tuple[int, int], ...]

    @property
    def coding_length(self):
        return sum(end - start for start, end in self.coding)


def read_models(path: Path) -> dict[str, Model]:
    models = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            starts = [int(v) for v in row["exonStarts"].strip(",").split(",")]
            ends = [int(v) for v in row["exonEnds"].strip(",").split(",")]
            if len(starts) != len(ends) or len(starts) != int(row["exonCount"]):
                raise ValueError("Exon count mismatch")
            exons = tuple(zip(starts, ends))
            if merge_blocks(exons) != exons:
                raise ValueError("Unsorted, overlapping or adjacent source exons")
            cds_start, cds_end = int(row["cdsStart"]), int(row["cdsEnd"])
            coding = tuple((max(s, cds_start), min(e, cds_end))
                           for s, e in exons if max(s, cds_start) < min(e, cds_end))
            if not coding or row["strand"] not in ("+", "-"):
                raise ValueError("Missing coding model")
            name = row["name"]
            if name in models:
                raise ValueError("Repeated source model name")
            model = Model(name, row["chrom"], row["strand"], exons, coding)
            encoded = re.fullmatch(
                r"PB\.\d+\.\d+:([^:|]+):([+-])\|(\d+)\|(\d+):(\d+):(\d+)\|([^|]+)\|([ACGT]{3})",
                name,
            )
            if name.startswith("PB."):
                if not encoded:
                    raise ValueError(f"Unrecognised PB source identity: {name}")
                chrom, strand, _, length, start, end, _, _ = encoded.groups()
                if chrom != model.chrom or strand != model.strand:
                    raise ValueError("Original ORF chromosome/strand mismatch")
                first = coding[0][0] if strand == "+" else coding[-1][1]
                last = coding[-1][1] if strand == "+" else coding[0][0]
                if (sum(e-s for s,e in exons) != int(length)
                    or boundary_offset(exons, strand, first) != int(start)-1
                    or boundary_offset(exons, strand, last) != int(end)-1):
                    raise ValueError(f"Transcript ORF coordinates disagree with genePred: {name}")
            models[name] = model
    return models


def check_gtf(path: Path, models: dict[str, Model]) -> dict:
    features = defaultdict(lambda: defaultdict(list))
    identities = defaultdict(set)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip("\n").split("\t")
            transcript = re.search(r'transcript_id "([^"]+)"', attributes)
            if transcript is None:
                raise ValueError("GTF row without transcript identity")
            name = transcript.group(1)
            # GTF is 1-based inclusive. Normalise ONCE to 0-based half-open.
            features[name][feature].append((int(start)-1, int(end)))
            identities[name].add((chrom, strand))
    if set(features) != set(models):
        raise ValueError("genePred and GTF model inventories differ")
    for name, model in models.items():
        f = features[name]
        if identities[name] != {(model.chrom, model.strand)}:
            raise ValueError("GTF chromosome/strand mismatch")
        if tuple(sorted(f["exon"])) != model.exons:
            raise ValueError(f"genePred/GTF exon mismatch: {name}")
        # genePred cdsEnd includes stop codon; this GTF lists it separately.
        if merge_blocks(f["CDS"] + f["stop_codon"]) != model.coding:
            raise ValueError(f"genePred/GTF coding mismatch: {name}")
        if sum(e-s for s,e in f["start_codon"]) != 3:
            raise ValueError("Incomplete source start codon")
    return {"models_checked": len(models), "all_exons_and_coding_blocks_equal": True,
            "GTF_conversion": "[start1,end1] -> [start1-1,end1)",
            "genePred_coding_includes_separate_GTF_stop_codon": True}


def join_feature(feature: dict, original_by_spacer: dict, models_by_id: dict,
                 enst_parents: dict) -> dict:
    spacer = reverse_complement(feature["sequence"].upper())
    original = original_by_spacer.get(spacer, [])
    if not original:
        raise ValueError(f"No original spacer for {feature['guide']}")
    def target_id(name):
        return enst_parents.get(name.split(".")[0], {}).get("Parent", identity(name))
    matched = [r for r in original if target_id(r["original_gene"]) == feature["target_gene_id"]]
    if len(matched) != 1:
        raise ValueError(f"Original target join is not unique: {feature['guide']}")
    row = matched[0]
    source_id = target_id(row["original_gene"])
    candidates = models_by_id.get(source_id, [])
    offsets = []
    for model in candidates:
        if model.chrom != row["chrom"]:
            continue
        offset = boundary_offset(model.coding, model.strand, int(row["cut0"]))
        offsets.append((model, offset))
    valid = [(model, offset) for model, offset in offsets if offset is not None]
    distances = {offset for _, offset in valid}
    base = {"guide": feature["guide"], "target_gene_id": feature["target_gene_id"],
            "original_gene": row["original_gene"], "spacer": spacer,
            "chrom": row["chrom"], "cut0": row["cut0"],
            "design_relative_position": row["relative"],
            "source_rows_with_this_spacer": [r["original_row"] for r in original],
            "source_aliases": [r["original_gene"] for r in original],
            "source_control_type": feature["source_control_type"],
            "RNA_QC_records": int(feature["single_guide_after_stage1_qc"])}
    if len(distances) != 1 or len(valid) != len(candidates):
        return dict(base, status="unresolved", distance=None,
                    reason="missing_or_outside_or_disagreeing_source_coding_model")
    distance = distances.pop()
    if distance != int(row["relative"]):
        raise ValueError(f"Recomputed distance disagrees with original design: {feature['guide']}: "
                         f"{distance} versus {row['relative']}")
    return dict(base, status="assessed", distance=distance,
                coding_lengths=sorted({m.coding_length for m,_ in valid}),
                annotation_models=[m.name for m,_ in valid],
                at_outer_ORF_boundary=any(distance in (0,m.coding_length) for m,_ in valid))


def make_plot(rows: list[dict], output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "svg.fonttype": "none", "axes.spines.top": False,
                         "axes.spines.right": False})
    distances = sorted(row["distance"] for row in rows if row["status"] == "assessed")
    if not distances:
        raise ValueError("No usable design distances; refusing an empty figure")
    counts = Counter(distances)
    x = sorted(counts)
    y, total = [], 0
    for distance in x:
        total += counts[distance]
        y.append(total)
    right = max(500, max(distances)) * 1.04
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    fig.subplots_adjust(left=.13, right=.96, bottom=.22, top=.83)
    ax.step([min(0,min(x)),*x,right], [0,*y,y[-1]], where="post",
            color="#367B97", linewidth=2.1)
    ax.axvline(0, color="#8F969B", linestyle="--", linewidth=1)
    ax.set_xlim(min(0,min(x))-right*.025, right)
    ax.set_ylim(0,len(rows)*1.04)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.set_xlabel("Predicted cut-site distance from CDS start (coding nt)", labelpad=10)
    ax.set_ylabel("Cumulative number of sgRNAs", labelpad=9)
    ax.grid(axis="y", color="#E6E8E9", linewidth=.6)
    ax.set_axisbelow(True)
    fig.suptitle("sgRNA design positions along the coding sequence", y=.965, fontsize=15)
    fig.text(.5,.885, f"Perturb-seq library | {len(rows):,} unique targeting sgRNAs; "
             f"{len(distances):,} with verified positions", ha="center", fontsize=11)
    fig.text(.13,.075,"Distance follows annotated coding exons; introns are excluded.", fontsize=10)
    fig.text(.13,.035,"All targeting designs included, regardless of RNA detection. These are predicted cut sites.",
             fontsize=9.5, color="#555B60")
    for suffix in ("png", "svg"):
        fig.savefig(output/f"design_cut_positions.{suffix}", dpi=220, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--original-sites", type=Path, required=True)
    parser.add_argument("--gene-table", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--enst-parents", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--code-commit", required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    models = read_models(args.gene_table)
    gtf_check = check_gtf(args.gtf, models)
    by_id = defaultdict(list)
    for model in models.values():
        by_id[identity(model.name)].append(model)
    source = json.loads(args.original_sites.read_text(encoding="utf-8"))
    by_spacer = defaultdict(list)
    for gene, sites in source["sites_by_original_Gene"].items():
        for site in sites:
            by_spacer[site["spacer"]].append(dict(site, original_gene=gene))
    parents = json.loads(args.enst_parents.read_text(encoding="utf-8"))["current_transcripts"]
    with args.features.open(encoding="utf-8-sig", newline="") as handle:
        features = list(csv.DictReader(handle))
    targeting = [row for row in features if row["target_gene_id"] != "Non-Targeting"]
    if len({row["sequence"] for row in targeting}) != len(targeting):
        raise ValueError("Repeated targeting sequence in expression feature inventory")
    if len(by_spacer) != len(targeting):
        raise ValueError("Original and observed feature design inventories differ")
    rows = [join_feature(row, by_spacer, by_id, parents) for row in targeting]
    assessed = [row for row in rows if row["status"] == "assessed"]
    distances = [row["distance"] for row in assessed]
    result = {
        "created_at": datetime.now(timezone.utc).isoformat(), "code_commit": args.code_commit,
        "source_files": {name:{"path":str(getattr(args,name)),"sha256":sha256(getattr(args,name))}
                         for name in ("features","original_sites","gene_table","gtf","enst_parents")},
        "method_sources": ["https://genome.ucsc.edu/FAQ/FAQformat.html#format9",
                           "https://genome.ucsc.edu/FAQ/FAQformat.html#format4"],
        "coordinate_contract": "genePred/cut0 0-based boundaries; GTF 1-based inclusive converted once; "
                               "sum coding exonic bases preceding cut along ORF strand, including source stop codon",
        "scope": "All unique targeting designs in the accepted feature inventory; no RNA/count/response filtering. "
                 "An intended cut-site audit, not an assay of actual indels or functional KO.",
        "gtf_gene_table_check": gtf_check,
        "original_rows": sum(map(len,source["sites_by_original_Gene"].values())),
        "original_unique_targeting_sequences": len(by_spacer),
        "duplicate_source_sequences": sum(len(v)>1 for v in by_spacer.values()),
        "designed_targeting_guides": len(rows), "assessed_guides": len(assessed),
        "unresolved_guides": len(rows)-len(assessed),
        "original_relative_positions_equal_for_all_assessed": True,
        "guides_with_zero_RNA_QC_records_retained": sum(row["RNA_QC_records"]==0 for row in rows),
        "summary": {"minimum_nt":min(distances),"maximum_nt":max(distances),
                    "median_nt":statistics.median(distances),
                    "within_nt":{str(n):sum(0<=d<=n for d in distances) for n in (50,100,200,300,500)},
                    "outer_ORF_boundary_guides":[row["guide"] for row in assessed if row["at_outer_ORF_boundary"]]},
        "records": rows,
    }
    (args.output/"design_cut_positions.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    make_plot(rows,args.output)
    print(json.dumps({k:v for k,v in result.items() if k not in ("records","source_files")},ensure_ascii=False))


if __name__ == "__main__":
    main()
