"""Render two meeting figures from the accepted RNA sequence and expression results."""
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


INK = "#23323B"
GRAY = "#D9E0E4"
BLUE = "#3C748E"
RED = "#AD5546"
PAGE = (11.2, 6.3)


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def variant_identity(row):
    return tuple(row[key] for key in ("guide", "chrom", "pos0", "ref", "alt"))


def sequence_counts(sequence_path: Path, coding_path: Path):
    source = read_json(sequence_path)["company"]
    coding = read_json(coding_path)
    guides = {row["guide"] for row in source["guides"] if row["primary_197"]}
    variants = [row for row in source["variants"] if row["primary_197"]]
    indels = [row for row in variants if row["kind"] in ("INS", "DEL")]
    assert len(guides) == coding["guide_count"]
    assert {variant_identity(row) for row in indels} == {
        variant_identity(row) for row in coding["events"]
    }
    assert len(indels) == len(coding["events"])
    assert all(row["guide"] in guides and row["guide_alt"] >= 1 for row in variants)
    support = Counter((row["kind"], row["guide_alt"] >= 2) for row in variants)
    effects = Counter((row["kind"], row["category"]) for row in coding["events"])
    assert all(category in ("Frameshift", "In-frame", "Unresolved") for _, category in effects)
    assert all(
        row["reason"] == "deletion_crosses_coding_or_splice_boundary"
        for row in coding["events"] if row["category"] == "Unresolved"
    )
    return support, effects, len(guides), coding["target_count"]


def expression_counts(root: Path):
    cases = read_csv(root / "case_summary.csv")
    results = {}
    for cohort in ("company", "chip01"):
        targets = {row["target_gene_id"] for row in cases if row["dataset"] == cohort}
        guides = read_json(root / f"{cohort}_profile_metadata.json")
        by_guide = {row["guide"]: row["target_gene_id"] for row in guides}
        assert len(guides) == len(by_guide)
        assert all(row["n_positive"] > 0 and row["n_reference"] > 0 for row in guides)
        assert set(by_guide.values()) <= targets
        by_target = Counter(by_guide.values())
        bins = [sum(min(by_target[target], 2) == n for target in targets) for n in range(3)]
        pairs = read_csv(root / f"{cohort}_guide_pairs.csv")
        expected = set(combinations(sorted(by_guide), 2))
        actual = {tuple(sorted((row["guide_left"], row["guide_right"]))) for row in pairs}
        assert len(pairs) == len(actual) and actual == expected
        same, different = [], []
        for row in pairs:
            same_target = by_guide[row["guide_left"]] == by_guide[row["guide_right"]]
            assert same_target == (row["same_target"].lower() == "true")
            value = float(row["correlation_raw"])
            assert np.isfinite(value) and -1 <= value <= 1
            (same if same_target else different).append(value)
        results[cohort] = (bins, same, different)
    return results


def title(fig, heading, subheading):
    fig.text(.065, .934, heading, fontsize=21, fontweight="bold", color=INK, va="top")
    fig.text(.065, .855, subheading, fontsize=11, color=INK, va="top", linespacing=1.45)


def panel_title(ax, label, pad=12):
    ax.set_title(label, loc="left", fontsize=12.5, fontweight="bold", pad=pad)


def count_bars(ax, labels, values, colors, legend_labels, *, show_totals=False):
    for y, row in zip(range(len(labels) - 1, -1, -1), values):
        left = 0
        for n, color in zip(row, colors):
            if n:
                ax.barh(y, n, left=left, height=.43, color=color, zorder=3)
                above = n < 7
                ax.text(left + n / 2, y + (.26 if above else 0), str(n), ha="center",
                        va="bottom" if above else "center", fontsize=10, zorder=4,
                        color=INK if above or color in (GRAY, "#A9C3D1") else "white")
            left += n
        if show_totals:
            ax.text(left + 2, y, str(left), va="center", fontsize=11, color=INK)
    ax.set_yticks(range(len(labels)), labels[::-1])
    ax.set_ylim(-.57, len(labels) - .37)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color="#E8ECEF", lw=.7)
    ax.legend(handles=[Patch(facecolor=color, label=label) for color, label in zip(colors, legend_labels)],
              loc="upper left", bbox_to_anchor=(0, -.27), ncol=len(colors),
              fontsize=9, frameon=False, borderaxespad=0, columnspacing=1.0, handlelength=1.1)


def sequence_figure(support, effects, guides, targets):
    fig = plt.figure(figsize=PAGE)
    title(fig, "Sequence changes observed in target RNA",
          f"{guides:,} sgRNAs across {targets} targets with detectable target RNA and reads in the target region")
    left, right = fig.add_gridspec(1, 2, left=.149, right=.963, bottom=.27, top=.731, wspace=.30).subplots()
    panel_title(left, "A   Sequence changes and RNA support")
    types = ("SNV", "INS", "DEL")
    count_bars(left, ["Base substitution", "Insertion", "Deletion"],
               [[support[kind, False], support[kind, True]] for kind in types],
               [GRAY, BLUE], ["1 UMI family", "≥2 UMI families"], show_totals=True)
    left.set_xlabel("Number of sequence changes", labelpad=8)
    panel_title(right, f"B   Predicted effect of {sum(effects.values())} indels")
    categories = ("Frameshift", "In-frame", "Unresolved")
    names = ("Frameshift", "In-frame", "Across ORF end")
    count_bars(right, ["Insertion", "Deletion"],
               [[effects[kind, category] for category in categories] for kind in ("INS", "DEL")],
               [RED, BLUE, GRAY],
               [f"{name} ({sum(effects[kind, category] for kind in ('INS', 'DEL'))})"
                for name, category in zip(names, categories)], show_totals=True)
    right.set_xlabel("Number of insertions and deletions", labelpad=8)
    limit = max(sum(support[kind, supported] for supported in (False, True)) for kind in types)
    for ax in (left, right):
        ax.set_xlim(0, limit * 1.14)
        ax.set_xticks(np.arange(0, limit + 1, 25))
    fig.text(.065, .054, "Each change is counted once per sgRNA. Reading-frame effects are predicted from the annotated ORF.",
             fontsize=9, color=INK)
    return fig


def expression_figure(results):
    fig = plt.figure(figsize=PAGE)
    title(fig, "Similarity of RNA expression changes across sgRNAs",
          "Within each sgRNA, compare RNA expression in records with predicted frameshift RNA\n"
          "and records with reference RNA at the corresponding site.")
    left, right = fig.add_gridspec(1, 2, width_ratios=(.87, 1.56), left=.127, right=.973,
                                  bottom=.26, top=.682, wspace=.47).subplots()
    panel_title(left, "A   Comparable sgRNAs per target", pad=24)
    cohorts = ("company", "chip01")
    count_bars(left, ["Company\nmerged", "Chip 1"], [results[cohort][0] for cohort in cohorts],
               [GRAY, "#A9C3D1", BLUE], ["0 sgRNAs", "1 sgRNA", "≥2 sgRNAs"])
    targets = max(sum(results[cohort][0]) for cohort in cohorts)
    left.set_xlim(0, targets * 1.035)
    left.set_xticks(np.arange(0, targets + 1, 25))
    left.set_xlabel("Number of targets", labelpad=8)

    panel_title(right, "B   Similarity of RNA changes", pad=24)
    positions = [0, 1.18, 2.8, 3.98]
    values = [data for cohort in cohorts for data in results[cohort][1:]]
    colors = [RED, "#A7B2B9", RED, "#A7B2B9"]
    boxes = right.boxplot(values, positions=positions, widths=.39, showfliers=False,
                          patch_artist=True, medianprops={"color": INK, "linewidth": 1.3},
                          whiskerprops={"color": "#667680"}, capprops={"color": "#667680"},
                          boxprops={"linewidth": .9}, zorder=3)
    for i, (box, pos, points, color) in enumerate(zip(boxes["boxes"], positions, values, colors)):
        box.set(facecolor="white", edgecolor=color, alpha=.8)
        y = np.sort(points)
        spread = (np.mod(np.arange(len(y)) * .61803398875, 1) - .5) * .38
        right.scatter(pos + spread, y, s=18 if i % 2 == 0 else 5,
                      c=color, alpha=.72 if i % 2 == 0 else .28, linewidths=0,
                      rasterized=False, zorder=4)
    right.set_xticks(positions, [f"{('Same target' if i % 2 == 0 else 'Different targets')}\n{len(data):,} pairs"
                                 for i, data in enumerate(values)], fontsize=9.5)
    for center, label in [(0.59, "Company merged"), (3.39, "Chip 1")]:
        right.text(center, 1.02, label, ha="center", va="bottom", fontsize=10,
                   fontweight="bold", transform=right.get_xaxis_transform())
    bound = max(.4, np.ceil(max(abs(r) for data in values for r in data) * 10) / 10)
    right.set_ylim(-bound, bound)
    right.set_yticks(np.linspace(-bound, bound, 5))
    right.set_ylabel("Pearson r of RNA changes", labelpad=8)
    right.axhline(0, color="#8B979E", linewidth=.9, zorder=1)
    right.grid(axis="y", color="#E8ECEF", lw=.7, zorder=0)
    right.set_xlim(-.62, 4.60)
    fig.text(.065, .092, f"Each source has {targets} targets. A sgRNA can be compared when both RNA groups are present.",
             fontsize=9, color=INK)
    fig.text(.065, .054, "Each point is one sgRNA pair. All eligible pairs are shown.", fontsize=9, color=INK)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sequences", "coding", "rna-results", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                         "svg.fonttype": "path"})
    sequence = sequence_counts(args.sequences, args.coding)
    expression = expression_counts(args.rna_results)
    figures = {"sequence_changes": sequence_figure(*sequence), "RNA_response_agreement": expression_figure(expression)}
    for name, fig in figures.items():
        for suffix in ("png", "svg"):
            fig.savefig(args.output / f"{name}.{suffix}", dpi=300, facecolor="white")
        plt.close(fig)
    inputs = [args.sequences, args.coding, args.rna_results / "case_summary.csv"]
    for cohort in ("company", "chip01"):
        inputs.extend(args.rna_results / f"{cohort}_{suffix}" for suffix in ("profile_metadata.json", "guide_pairs.csv"))
    receipt = {"code_commit": args.commit,
               "inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
               "figures_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.iterdir()},
               "target_counts_by_0_1_ge2_comparable_guides": {k: v[0] for k, v in expression.items()},
               "plotted_pair_counts": {k: {"same": len(v[1]), "different": len(v[2])} for k, v in expression.items()}}
    (args.output / "render_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ("target_counts_by_0_1_ge2_comparable_guides", "plotted_pair_counts")}, indent=2))


if __name__ == "__main__":
    main()
