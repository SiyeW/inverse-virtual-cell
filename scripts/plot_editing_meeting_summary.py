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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from render_share_ready_figures import COLORS, LABELS, assert_drawn_points, configure, save


GRAY = "#BAC2C7"
BLUE = "#417893"
FRAME = "#A47C55"


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
    return support, effects


def guide_sequence_summary(sequence_path: Path, coding_path: Path):
    source = read_json(sequence_path)["company"]
    coding = read_json(coding_path)
    guides = {row["guide"] for row in source["guides"] if row["primary_197"]}
    assert len(guides) == coding["guide_count"] == 197
    variants = [row for row in source["variants"] if row["primary_197"]]
    members = {"Any change": {row["guide"] for row in variants}}
    for label, kind in (("Substitution", "SNV"), ("Insertion", "INS"), ("Deletion", "DEL")):
        members[label] = {row["guide"] for row in variants if row["kind"] == kind}
    members["Any indel"] = {row["guide"] for row in coding["events"]}
    for label, category in (("Frameshift", "Frameshift"), ("In-frame", "In-frame"),
                            ("Across ORF end", "Unresolved")):
        members[label] = {row["guide"] for row in coding["events"] if row["category"] == category}
    assert members["Any indel"] == members["Insertion"] | members["Deletion"]
    assert all(group <= guides for group in members.values())
    assert members["Any change"] == members["Substitution"] | members["Any indel"]
    return {"denominator": len(guides), "counts": {key: len(group) for key, group in members.items()},
            "guide_membership": {key: sorted(group) for key, group in members.items()}}


def guide_sequence_figure(summary, *, percent=True):
    denominator = summary["denominator"]
    counts = summary["counts"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.98, bottom=.18, top=.9, wspace=.19)
    definitions = [(["Any change", "Substitution", "Insertion", "Deletion"],
                    [BLUE] * 4, "RNA sequence changes"),
                   (["Frameshift", "In-frame", "Across ORF end"],
                    [FRAME, BLUE, "white"], "Predicted indel effects")]
    for ax, (labels, colors, title) in zip(axes, definitions):
        values = [counts[label] / denominator * 100 if percent else counts[label] for label in labels]
        bars = ax.bar(range(len(labels)), values, width=.56, color=colors, linewidth=.6,
                      edgecolor=["#66747C" if color == "white" else "white" for color in colors])
        for bar, label, value in zip(bars, labels, values):
            if label == "Across ORF end":
                bar.set_hatch("///")
            text = f"{value:.1f}%\n{counts[label]}/{denominator}" if percent else str(counts[label])
            ax.text(bar.get_x() + bar.get_width() / 2, value + (1.2 if percent else 2),
                    text, ha="center", va="bottom", fontsize=9, linespacing=1.25)
        assert np.allclose([bar.get_height() for bar in bars], values, rtol=0, atol=1e-12)
        ax.set_xticks(range(len(labels)), [label.replace("Across ORF end", "Across\nORF end") for label in labels])
        ax.set_title(title)
        ax.set_xlim(-.6, len(labels) - .4)
    axes[0].set_ylabel(f"sgRNAs (% of {denominator})" if percent else "Number of sgRNAs")
    axes[0].set_ylim(0, 70 if percent else 140)
    axes[0].set_yticks(np.arange(0, 71, 10) if percent else np.arange(0, 141, 20))
    return fig


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


def count_bars(ax, labels, values, styles, legend_labels, *, show_totals=False, percent=False):
    expected_heights = []
    for x, row in enumerate(values):
        total = sum(row)
        assert total > 0
        bottom = 0
        for n, style in zip(row, styles):
            assert isinstance(n, int) and n >= 0
            if n:
                height = n / total * 100 if percent else n
                ax.bar(x, height, bottom=bottom, width=.56, linewidth=.55, **style)
                expected_heights.append(height)
                small = height < 7
                light = style["facecolor"] in (GRAY, "white")
                label = f"{height:.1f}%" if percent else str(n)
                ax.text(x, bottom + (height + 1.9 if small else height / 2), label, ha="center",
                        va="bottom" if small else "center", fontsize=9,
                        color="#333333" if small or light else "white")
                bottom += height
        assert np.isclose(bottom, 100 if percent else total)
        if show_totals:
            ax.text(x, 110 if percent else bottom + 2.1,
                    f"n = {total}" if percent else str(total), ha="center", va="bottom", fontsize=9)
    assert np.allclose([patch.get_height() for patch in ax.patches], expected_heights, rtol=0, atol=1e-12)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_xlim(-.54, len(labels) - .46)
    ax.legend(handles=[Patch(label=label, **style) for style, label in zip(styles, legend_labels)],
              loc="upper center", bbox_to_anchor=(.5, -.14), ncol=len(styles),
              columnspacing=1.1, handlelength=1.4, handletextpad=.45)


def sequence_figure(support, effects, *, percent=False):
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 4.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.98, bottom=.22, top=.91, wspace=.19)
    left.set_title("Sequence changes in target RNA")
    gray = {"facecolor": GRAY, "edgecolor": "white"}
    blue = {"facecolor": BLUE, "edgecolor": "white"}
    types = ("SNV", "INS", "DEL")
    count_bars(left, ["Substitution", "Insertion", "Deletion"],
               [[support[kind, False], support[kind, True]] for kind in types],
               [gray, blue], ["1 UMI family", "≥2 UMI families"], show_totals=True, percent=percent)
    left.set_ylabel("Proportion of sequence changes (%)" if percent else "Number of RNA sequence changes")
    right.set_title("Predicted effect of insertions and deletions")
    categories = ("Frameshift", "In-frame", "Unresolved")
    count_bars(right, ["Insertion", "Deletion"],
               [[effects[kind, category] for category in categories] for kind in ("INS", "DEL")],
               [{"facecolor": FRAME, "edgecolor": "white"}, blue,
                {"facecolor": "white", "edgecolor": "#66747C", "hatch": "///"}],
               ["Frameshift", "In-frame", "Across ORF end"], show_totals=percent, percent=percent)
    limit = max(sum(support[kind, supported] for supported in (False, True)) for kind in types)
    left.set_ylim(0, 117 if percent else limit * 1.14)
    left.set_yticks(np.arange(0, 101, 25) if percent else np.arange(0, limit + 1, 25))
    return fig


def expression_figure(results):
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    fig.subplots_adjust(left=.14, right=.98, bottom=.13, top=.96)
    expected = []
    for j, (cohort, color) in enumerate(zip(("chip01", "company"), COLORS)):
        for i, points in enumerate(results[cohort][1:]):
            position = i + (j - .5) * .24
            box = ax.boxplot([points], positions=[position], widths=.19, patch_artist=True,
                             manage_ticks=False, showfliers=False,
                             medianprops={"color": color, "linewidth": 1.25},
                             whiskerprops={"color": color, "linewidth": .75},
                             capprops={"color": color, "linewidth": .75})
            box["boxes"][0].set(facecolor=matplotlib.colors.to_rgba(color, .12),
                                edgecolor=color, linewidth=.75)
            x = position + np.random.default_rng(90210 + j * 10 + i).uniform(-.065, .065, len(points))
            ax.scatter(x, points, s=5.5, alpha=.35, color=color, edgecolors="none", zorder=2)
            expected.extend(zip(x, points))
    assert_drawn_points(ax, expected)
    handles = [Line2D([], [], marker="o", markersize=4, linestyle="none", color=color, label=label)
               for color, label in zip(COLORS[:2], LABELS[:2])]
    ax.legend(handles=handles, loc="upper right", ncol=2, columnspacing=1.3, handletextpad=.4)
    ax.set_xticks([0, 1], ["Same target", "Different targets"])
    bound = max(.4, np.ceil(max(abs(r) for _, r in expected) * 10) / 10)
    ax.set(xlim=(-.48, 1.48), ylim=(-bound, bound), ylabel="Correlation of RNA changes (Pearson r)")
    ax.set_yticks(np.linspace(-bound, bound, 5))
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#ededed", linewidth=.45)
    ax.axhline(0, color="#999999", linewidth=.65, zorder=1)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sequences", "coding", "rna-results", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--sequence-scale", choices=("count", "percent"), default="count")
    parser.add_argument("--sequence-unit", choices=("change", "guide"), default="change")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    font = Path("C:/Windows/Fonts/arial.ttf")
    configure(font if font.is_file() else None)
    sequence = sequence_counts(args.sequences, args.coding)
    guide_summary = guide_sequence_summary(args.sequences, args.coding)
    expression = expression_counts(args.rna_results)
    sequence_plot = (guide_sequence_figure(guide_summary, percent=args.sequence_scale == "percent")
                     if args.sequence_unit == "guide"
                     else sequence_figure(*sequence, percent=args.sequence_scale == "percent"))
    figures = {"sequence_changes": sequence_plot,
               "RNA_response_agreement": expression_figure(expression)}
    for name, fig in figures.items():
        save(fig, args.output, name)
    inputs = [args.sequences, args.coding, args.rna_results / "case_summary.csv"]
    for cohort in ("company", "chip01"):
        inputs.extend(args.rna_results / f"{cohort}_{suffix}" for suffix in ("profile_metadata.json", "guide_pairs.csv"))
    receipt = {"code_commit": args.commit,
               "sequence_scale": args.sequence_scale,
               "sequence_unit": args.sequence_unit,
               "guide_sequence_summary": guide_summary,
               "inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
               "figures_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.iterdir()},
               "sequence_support_counts": {f"{kind}_{'ge2' if supported else 'one'}": n
                                           for (kind, supported), n in sequence[0].items()},
               "indel_consequence_counts": {f"{kind}_{category}": n
                                             for (kind, category), n in sequence[1].items()},
               "all_scatter_values_and_bar_heights_verified": True,
               "target_counts_by_0_1_ge2_comparable_guides": {k: v[0] for k, v in expression.items()},
               "plotted_pair_counts": {k: {"same": len(v[1]), "different": len(v[2])} for k, v in expression.items()}}
    (args.output / "render_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ("target_counts_by_0_1_ge2_comparable_guides", "plotted_pair_counts")}, indent=2))


if __name__ == "__main__":
    main()
