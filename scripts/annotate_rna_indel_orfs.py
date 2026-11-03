"""Predict individual RNA-observed indel effects against explicit ORF references.

This does not infer DNA genotypes, complete edited transcripts, or functional KO.
Every compatible reference is evaluated; disagreement remains unresolved.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path


CATEGORIES = ("Frameshift", "In-frame", "Outside ORF", "Unresolved")
COLORS = ("#BE7953", "#447E9C", "#E1E4E6", "#8C959C")


def rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@dataclass(frozen=True)
class ORF:
    chrom: str
    strand: str
    blocks: tuple[tuple[int, int], ...]
    sequence: str

    def __post_init__(self):
        assert self.strand in ("+", "-") and self.blocks
        assert all(0 <= start < end for start, end in self.blocks)
        assert all(a[1] <= b[0] for a, b in zip(self.blocks, self.blocks[1:]))
        assert sum(end - start for start, end in self.blocks) == len(self.sequence)
        assert len(self.sequence) % 3 == 0
        assert set(self.sequence) <= set("ACGT")

    def base_offset(self, position: int) -> int | None:
        prefix = 0
        ordered = self.blocks if self.strand == "+" else reversed(self.blocks)
        for start, end in ordered:
            if start <= position < end:
                return prefix + (position - start if self.strand == "+" else end - 1 - position)
            prefix += end - start
        return None

    def genomic_base(self, position: int) -> str | None:
        offset = self.base_offset(position)
        if offset is None:
            return None
        base = self.sequence[offset]
        return base if self.strand == "+" else rc(base)


def minimal_change(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Remove unchanged VCF anchors without shifting across repeats."""
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt = ref[prefix:], alt[prefix:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return position + prefix, ref, alt


def classify_one(event: dict, model: ORF) -> tuple[str, str, bool]:
    if event["chrom"] != model.chrom:
        return "Unresolved", "reference_chromosome_mismatch", False
    ref, alt = event["ref"].upper(), event["alt"].upper()
    if not ref or not alt or not set(ref + alt) <= set("ACGT"):
        return "Unresolved", "nonliteral_or_ambiguous_allele", False
    start, deleted, inserted = minimal_change(int(event["pos0"]), ref, alt)
    delta = len(inserted) - len(deleted)
    if delta == 0:
        return "Unresolved", "not_a_length_changing_indel", False
    if int(event["length_change"]) != delta:
        raise ValueError("Indel length disagrees with the accepted source event")
    if deleted:
        offsets = [model.base_offset(p) for p in range(start, start + len(deleted))]
        if all(offset is None for offset in offsets):
            return "Outside ORF", "deleted_bases_outside_reference_coding_blocks", False
        if any(offset is None for offset in offsets):
            return "Unresolved", "deletion_crosses_coding_or_splice_boundary", False
        observed_ref = "".join(model.genomic_base(p) for p in range(start, start + len(deleted)))
        if observed_ref != deleted:
            return "Unresolved", "source_REF_does_not_match_ORF_reference", False
        edge = any(offset < 3 or offset >= len(model.sequence) - 3 for offset in offsets)
    else:
        before, after = model.base_offset(start - 1), model.base_offset(start)
        if before is None or after is None:
            if start in (model.blocks[0][0], model.blocks[-1][1]):
                return "Outside ORF", "insertion_at_outer_reference_ORF_boundary", False
            if any(start in block for block in model.blocks):
                return "Unresolved", "insertion_at_internal_splice_boundary", False
            return "Outside ORF", "insertion_outside_reference_coding_blocks", False
        if abs(before - after) != 1:
            return "Unresolved", "insertion_flanks_not_adjacent_in_coding_sequence", False
        # Accept either VCF anchoring direction; validate the unchanged REF bases
        # that are represented by this coding reference.
        for i, base in enumerate(ref):
            expected = model.genomic_base(int(event["pos0"]) + i)
            if expected is not None and expected != base:
                return "Unresolved", "source_REF_does_not_match_ORF_reference", False
        edge = min(before, after) < 3 or max(before, after) >= len(model.sequence) - 3
    return ("In-frame" if delta % 3 == 0 else "Frameshift"), "coding_length_change", edge


def classify(event: dict, models: list[ORF]) -> dict:
    if not models:
        return {"category": "Unresolved", "reason": "target_ORF_reference_not_recovered", "start_or_stop_overlap": False}
    calls = [classify_one(event, model) for model in models]
    categories = {call[0] for call in calls}
    if len(categories) != 1:
        category, reason = "Unresolved", "compatible_ORF_references_disagree"
    else:
        category = calls[0][0]
        reason = ";".join(sorted({call[1] for call in calls}))
    return {"category": category, "reason": reason, "start_or_stop_overlap": any(call[2] for call in calls)}


def references(packet: dict) -> dict[str, list[ORF]]:
    if not packet["all_197_original_23nt_targets_match_GRCh38"]:
        raise ValueError("Original guide identity/reference validation has not passed")
    grouped = defaultdict(dict)
    for row in packet["models"]:
        if not row["reference_sequence_consistent_with_design"]:
            continue
        sequence = row["reference_ORF_sequence"]
        assert hashlib.sha256(sequence.encode()).hexdigest() == row["reference_sequence_sha256"]
        model = ORF(row["chrom"], row["strand"], tuple(map(tuple, row["validated_blocks"])), sequence)
        grouped[row["target_gene_id"]][model] = model
    return {gene: list(models) for gene, models in grouped.items()}


def analyze(sequence_packet: dict, reference_packet: dict, expected_guides: int, expected_indels: int) -> dict:
    cohort = [r for r in sequence_packet["company"]["guides"] if r["primary_197"]]
    assert len(cohort) == len({r["guide"] for r in cohort}) == expected_guides
    guide_genes = {r["guide"]: r["target_gene_id"] for r in cohort}
    models = references(reference_packet)
    events = [r for r in sequence_packet["company"]["variants"] if r["primary_197"] and r["kind"] in ("INS", "DEL")]
    assert len(events) == expected_indels
    identities = {(r["guide"], r["chrom"], r["pos0"], r["ref"], r["alt"]) for r in events}
    assert len(identities) == len(events)
    results = []
    for event in events:
        gene = guide_genes[event["guide"]]
        assert event["target_gene_id"] == gene
        keep = {key: event[key] for key in (
            "guide", "target_gene_id", "chrom", "pos0", "ref", "alt", "kind", "length_change",
            "guide_alt", "guide_ref", "NTC_alt", "NTC_ref",
        )}
        keep.update(classify(event, models.get(gene, [])))
        keep["compatible_reference_count"] = len(models.get(gene, []))
        results.append(keep)
    return {
        "guide_count": len(cohort), "target_count": len(set(guide_genes.values())),
        "indel_observations": len(results),
        "counting_unit": "Each distinct normalized indel counted once per sgRNA; not DNA alleles, independent cells or KO events",
        "reference_targets": len(models),
        "guides_with_reference": sum(gene in models for gene in guide_genes.values()),
        "counts": dict(Counter(r["category"] for r in results)),
        "reasons": dict(Counter(r["reason"] for r in results)),
        "initiation_or_termination_overlap_observations": sum(r["start_or_stop_overlap"] for r in results),
        "events": results,
    }


def plot(result: dict, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(9.0, 4.5))
    for y, kind in [(1, "INS"), (0, "DEL")]:
        counts = Counter(r["category"] for r in result["events"] if r["kind"] == kind)
        left = 0
        for category, color in zip(CATEGORIES, COLORS):
            value = counts[category]
            ax.barh(y, value, left=left, height=0.44, color=color, edgecolor="white", linewidth=0.7)
            if value >= 3:
                ax.text(left + value / 2, y, str(value), ha="center", va="center", color="black" if category == "Outside ORF" else "white")
            left += value
        ax.text(left + 1.2, y, f"n = {left}", ha="left", va="center")
    ax.set_yticks([1, 0], ["Insertions", "Deletions"])
    ax.set_ylim(-0.65, 1.65)
    ax.set_xlim(0, max(sum(r["kind"] == kind for r in result["events"]) for kind in ["INS", "DEL"]) * 1.15)
    ax.set_xlabel("Indel observations across sgRNAs")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=10)
    handles = [Patch(facecolor=c, label=f"{name} (n = {result['counts'].get(name, 0)})") for name, c in zip(CATEGORIES, COLORS)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.52, 0.86), ncol=2, frameon=False)
    fig.suptitle("Predicted coding effects of RNA-observed indels", y=0.98, fontsize=14)
    fig.text(0.13, 0.07, "Each distinct indel is counted once per sgRNA. All 197 sgRNAs are retained in the input cohort.", fontsize=9)
    fig.text(0.13, 0.025, "Unresolved: no verified reference, a boundary change, or reference disagreement. This is not KO efficiency.", fontsize=9)
    fig.subplots_adjust(left=0.13, right=0.91, top=0.66, bottom=0.23)
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--expected-guides", type=int, default=197)
    parser.add_argument("--expected-indels", type=int, default=145)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise FileExistsError("Use a new output directory; preserve previous accepted results")
    source = json.loads(args.sequences.read_text(encoding="utf-8-sig"))
    ref = json.loads(args.references.read_text(encoding="utf-8-sig"))
    result = analyze(source, ref, args.expected_guides, args.expected_indels)
    result["analysis_commit"] = args.commit
    result["input_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.sequences, args.references]}
    (args.output / "coding_consequences.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    plot(result, args.output / "coding_consequences")
    print(json.dumps({k: v for k, v in result.items() if k != "events"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
