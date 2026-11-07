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


def count_bars(ax, labels, values, styles, legend_labels, *, show_totals=False):
    expected_heights = []
    for x, row in enumerate(values):
        bottom = 0
        for n, style in zip(row, styles):
            assert isinstance(n, int) and n >= 0
            if n:
                ax.bar(x, n, bottom=bottom, width=.56, linewidth=.55, **style)
                expected_heights.append(n)
                small = n < 7
                light = style["facecolor"] in (GRAY, "white")
                ax.text(x, bottom + (n + 1.9 if small else n / 2), str(n), ha="center",
                        va="bottom" if small else "center", fontsize=9,
                        color="#333333" if small or light else "white")
            bottom += n
        if show_totals:
            ax.text(x, bottom + 2.1, str(bottom), ha="center", va="bottom", fontsize=9)
    assert [patch.get_height() for patch in ax.patches] == expected_heights
    ax.set_xticks(range(len(labels)), labels)
    ax.set_xlim(-.54, len(labels) - .46)
    ax.legend(handles=[Patch(label=label, **style) for style, label in zip(styles, legend_labels)],
              loc="upper center", bbox_to_anchor=(.5, -.14), ncol=len(styles),
              columnspacing=1.1, handlelength=1.4, handletextpad=.45)


def sequence_figure(support, effects):
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 4.8), sharey=True)
    fig.subplots_adjust(left=.09, right=.98, bottom=.22, top=.91, wspace=.19)
    left.set_title("Sequence changes in target RNA")
    gray = {"facecolor": GRAY, "edgecolor": "white"}
    blue = {"facecolor": BLUE, "edgecolor": "white"}
    types = ("SNV", "INS", "DEL")
    count_bars(left, ["Substitution", "Insertion", "Deletion"],
               [[support[kind, False], support[kind, True]] for kind in types],
               [gray, blue], ["1 UMI family", "≥2 UMI families"], show_totals=True)
    left.set_ylabel("Number of RNA sequence changes")
    right.set_title("Predicted effect of insertions and deletions")
    categories = ("Frameshift", "In-frame", "Unresolved")
    count_bars(right, ["Insertion", "Deletion"],
               [[effects[kind, category] for category in categories] for kind in ("INS", "DEL")],
               [{"facecolor": FRAME, "edgecolor": "white"}, blue,
                {"facecolor": "white", "edgecolor": "#66747C", "hatch": "///"}],
               ["Frameshift", "In-frame", "Across ORF end"])
    limit = max(sum(support[kind, supported] for supported in (False, True)) for kind in types)
    left.set_ylim(0, limit * 1.14)
    left.set_yticks(np.arange(0, limit + 1, 25))
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
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    font = Path("C:/Windows/Fonts/arial.ttf")
    configure(font if font.is_file() else None)
    sequence = sequence_counts(args.sequences, args.coding)
    expression = expression_counts(args.rna_results)
    figures = {"sequence_changes": sequence_figure(*sequence), "RNA_response_agreement": expression_figure(expression)}
    for name, fig in figures.items():
        save(fig, args.output, name)
    inputs = [args.sequences, args.coding, args.rna_results / "case_summary.csv"]
    for cohort in ("company", "chip01"):
        inputs.extend(args.rna_results / f"{cohort}_{suffix}" for suffix in ("profile_metadata.json", "guide_pairs.csv"))
    receipt = {"code_commit": args.commit,
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
