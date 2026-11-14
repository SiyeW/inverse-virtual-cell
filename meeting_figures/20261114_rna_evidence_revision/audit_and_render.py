"""Read saved small results, audit RNA-family units, and redraw labels only.

No BAM, expression matrix, model fit, DEG test, PCA fit or score is run.
All writes are confined to this new directory; source hashes are checked again.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import ast
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import posixpath
import subprocess
import sys
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.text import Text
import numpy as np
from PIL import Image, ImageOps, ImageDraw

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROJECT = REPO.parent
OPS = PROJECT / "operations"
SEQUENCE = OPS / "20261103_editing_cell_exploration/sequence_export_r2.json"
LABELS = OPS / "20261105_frameshift_RNA_states/coding_rna_record_labels.json"
CATALOGUE = OPS / "20261030_guide_rna_difference/catalogue.csv"
WORKBOOK = Path("D:/SiyeWu/Documents/xwechat_files/wxid_yrlf5y07mu3b12_2e89/msg/file/2026-08/crop-seq-sgRNA文库-v2.xlsx")
CODING = OPS / "20261104_project_ORF_indels/consequences_r1/coding_consequences.json"
COORDS = OPS / "20261110_shared_PCA_preview/shared_pca_coordinates.npz"
PPT = Path("D:/BGM/261107/de novo组会261107.pptx")
BLUE, GREY, BROWN, RED = "#417893", "#BAC2C7", "#A47C55", "#D33E35"
INPUTS = {}
FIGURES = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def track(path):
    path = Path(path)
    INPUTS[str(path)] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    return path


def read_json(path):
    return json.loads(track(path).read_text(encoding="utf-8-sig"))


def write_json(name, value):
    (HERE / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def table(name, rows):
    rows = list(rows)
    assert rows
    with (HERE / name).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def module(name, path):
    track(path)
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def style():
    plt.rcParams.update({"font.family": "Arial", "font.size": 11, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "axes.unicode_minus": False, "figure.facecolor": "white", "legend.frameon": False})


def save(fig, stem, caption, sources):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for artist in fig.findobj(Text):
        if not artist.get_visible() or not artist.get_text() or (artist.axes is not None and not artist.axes.get_visible()):
            continue
        bounds = artist.get_window_extent(renderer)
        if bounds.width and bounds.height and (bounds.x0 < -2 or bounds.y0 < -2 or
                bounds.x1 > fig.bbox.width + 2 or bounds.y1 > fig.bbox.height + 2):
            outside.append(artist.get_text())
    assert not outside, (stem, outside)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    root = ET.parse(HERE / f"{stem}.svg").getroot()
    assert not any(x.tag.endswith("}image") for x in root.iter()), "Raster in editable SVG"
    with Image.open(HERE / f"{stem}.png") as im:
        size = im.size
    FIGURES.append({"stem": stem, "pixels": size, "caption": caption,
        "sources": [str(p) for p in sources], "png_sha256": sha(HERE / f"{stem}.png"),
        "svg_sha256": sha(HERE / f"{stem}.svg"), "outside_text": outside})
    print("SAVED", stem, size, flush=True)


def workbook_rows():
    """Only parse the design's XML cells; no workbook edits or recalculation."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(track(WORKBOOK)) as z:
        shared = ["".join(x.itertext()) for x in ET.fromstring(z.read("xl/sharedStrings.xml"))]
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet.findall(".//m:sheetData/m:row", ns)[1:]:
            values = {}
            for cell in row:
                node = cell.find("m:v", ns)
                text = node.text if node is not None else ""
                key = "".join(x for x in cell.get("r") if x.isalpha())
                values[key] = shared[int(text)] if cell.get("t") == "s" else text
            rows.append(values)
    return rows


def design_audit():
    rows = list(csv.DictReader(track(CATALOGUE).open(encoding="utf-8-sig")))
    groups = defaultdict(list)
    for row in rows:
        if row["source_control_type"] in ("targeting", "positive_control"):
            groups[row["target_gene_id"], row["target_gene_name"]].append(row)
    target_rows = [{"target_gene_id": key[0], "target_gene_name": key[1],
        "catalogue_distinct_guides": len(value),
        "guides_with_stage1_QC_records": sum(int(x["single_guide_after_stage1_qc"]) > 0 for x in value),
        "stage1_QC_records": sum(int(x["single_guide_after_stage1_qc"]) for x in value),
        "guide_names": ";".join(x["guide"] for x in value)} for key, value in sorted(groups.items())]
    table("catalogue_target_guide_counts.csv", target_rows)
    table("pilot_guide_counts.csv", [r for r in rows if r["target_gene_name"] in ("TYMSOS", "UBA1")])
    original = workbook_rows()
    control_names = {"nontarget control", "nontarget-2020Science", "intergenic-2020Science", "intergenic"}
    targets = [r for r in original if r["A"] not in control_names]
    counts = Counter(r["A"] for r in targets)
    table("original_design_target_ID_counts.csv", [{"original_target_ID": k, "design_rows": v} for k, v in sorted(counts.items())])
    by_seq = defaultdict(list)
    for row in original:
        by_seq[row["D"]].append(row)
    duplicate_rows = [{"spacer_sequence": s, "design_target_IDs": ";".join(r["A"] for r in rr),
                       "rows": len(rr)} for s, rr in by_seq.items() if len(rr) > 1]
    table("original_design_duplicate_sequences.csv", duplicate_rows)
    assert len(original) == 1201 and len(by_seq) == len(rows) == 1196
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    catalogue_seqs = {r["sequence"] for r in rows}
    assert all(r["D"] in catalogue_seqs or rc(r["D"]) in catalogue_seqs for r in original)
    assert len(duplicate_rows) == 5 and all(r["rows"] == 2 for r in duplicate_rows)
    result = {"original_design_rows": len(original), "original_unique_spacers": len(by_seq),
        "original_target_IDs": len(counts), "original_design_target_ID_guide_distribution": dict(sorted(Counter(counts.values()).items())),
        "catalogue_rows": len(rows), "catalogue_targets": len(groups),
        "catalogue_target_guide_distribution": dict(sorted(Counter(len(v) for v in groups.values()).items())),
        "postQC_observed_guide_distribution": dict(sorted(Counter(r["guides_with_stage1_QC_records"] for r in target_rows).items())),
        "duplicates": duplicate_rows, "pilot": [r for r in target_rows if r["target_gene_name"] in ("TYMSOS", "UBA1")]}
    write_json("design_audit.json", result)
    return result


def identity(v):
    return tuple(v[k] for k in ("guide", "chrom", "pos0", "ref", "alt"))


def evidence_audit(data, labels):
    events = [v for v in data["variants"] if v["target"] in ("TYMSOS", "UBA1") and v["kind"] in ("INS", "DEL")]
    canonical = {identity(v): v for v in labels["variants"]}
    supported, callables, event_rows, family_rows = [], [], [], []
    displays = {d["guide"]: d for d in data["display"]}
    for v in events:
        original = canonical[identity(v)]
        assert original["cell_support"] == v["cell_support"]
        for k in ("alt", "ref", "other", "ambiguous"):
            assert sum(c[k] for c in v["cell_support"]) == v["guide_" + k]
        counts = {k: v["guide_" + k] for k in ("alt", "ref", "other", "ambiguous")}
        positives = [c for c in v["cell_support"] if c["alt"] > 0]
        event_rows.append({"target": v["target"], "guide": v["guide"], "chrom": v["chrom"],
            "position_1based": v["pos0"] + 1, "REF_allele": v["ref"], "ALT_allele": v["alt"],
            "kind": v["kind"], "length_change": v["length_change"],
            **{k + "_UMI_families": n for k, n in counts.items()},
            "informative_UMI_families": counts["alt"] + counts["ref"] + counts["other"],
            "ALT_support_record_n": len(positives), "REF_support_record_n": sum(c["ref"] > 0 for c in v["cell_support"]),
            "scan_capped": v["scan_capped"]})
        for cell in v["cell_support"]:
            total = sum(cell[k] for k in ("alt", "ref", "other"))
            r = {"dataset": "company", "target": v["target"], "guide": v["guide"],
                "cell_barcode": cell["cell_barcode"], "chrom": v["chrom"], "position_1based": v["pos0"] + 1,
                "REF_allele": v["ref"], "ALT_allele": v["alt"],
                **{k + "_UMI_families": cell[k] for k in ("alt", "ref", "other", "ambiguous")},
                "informative_UMI_families": total,
                "ALT_fraction_observed": cell["alt"] / total if total else "",
                "ALT_and_REF_observed": bool(cell["alt"] and cell["ref"]),
                "scan_capped": v["scan_capped"],
                "coverage_interpretation": "Only 1-2 informative families; cannot exclude unobserved REF" if cell["alt"] else ""}
            callables.append(r)
            if cell["alt"]:
                supported.append(r)
    table("event_counts_and_PCA_n.csv", event_rows)
    table("all_event_site_record_counts.csv", callables)
    table("indel_supported_record_evidence_40.csv", supported)
    table("two_example_record_counts.csv", [r for r in callables if
        (r["guide"] == "TYMSOS-2" and r["position_1based"] == 658115) or
        (r["guide"] == "UBA1-2" and r["position_1based"] == 47198890)])
    wanted = {(r["guide"], r["cell_barcode"]) for r in supported}
    for d in displays.values():
        for m in d["molecules"]:
            if (d["guide"], m["cell_barcode"]) in wanted:
                fields = m["SAM"].split("\t")
                family_rows.append({"guide": d["guide"], "cell_barcode": m["cell_barcode"],
                    "corrected_UMI": m["UMI"], "GX": m["GX"],
                    "saved_deduplicated_alignment_observations": m["PCR_observations"],
                    "representative_start_1based": fields[3], "representative_CIGAR": fields[5],
                    "note": "Regional representative only; exact family-state counts come from cell_support"})
    table("support_record_regional_family_identifiers.csv", family_rows)
    duplicates = Counter((r["guide"], r["cell_barcode"], r["corrected_UMI"]) for r in family_rows)
    assert not any(n > 1 for n in duplicates.values())
    summaries = []
    for target, expected in (("TYMSOS", 9), ("UBA1", 31)):
        rr = [r for r in supported if r["target"] == target]
        assert len(rr) == len({(r["guide"], r["cell_barcode"]) for r in rr}) == expected
        summaries.append({"target": target, "indel_support_records": len(rr),
            "ALT_and_REF_observed_records": sum(r["ALT_and_REF_observed"] for r in rr),
            "only_ALT_observed_records": sum(r["ref_UMI_families"] == r["other_UMI_families"] == r["ambiguous_UMI_families"] == 0 for r in rr),
            "one_informative_family_records": sum(r["informative_UMI_families"] == 1 for r in rr),
            "two_informative_family_records": sum(r["informative_UMI_families"] == 2 for r in rr),
            "total_ALT_families": sum(r["alt_UMI_families"] for r in rr),
            "zero_callable_records_among_ALT_supported": 0,
            "records_insufficient_to_exclude_REF": len(rr)})
    table("record_coverage_summary.csv", summaries)
    dd = displays["TYMSOS-2"]
    deletions = sorted([v for v in events if v["guide"] == "TYMSOS-2" and v["kind"] == "DEL"], key=lambda v: v["pos0"])
    alleles = []
    for v in deletions:
        offset = v["pos0"] - dd["start0"]
        assert dd["reference"][offset:offset+len(v["ref"])] == v["ref"]
        alleles.append({"position_1based": v["pos0"] + 1, "REF": v["ref"], "ALT": v["alt"],
            "mutated_common_window": dd["reference"][:offset] + v["alt"] + dd["reference"][offset+len(v["ref"]):]})
    assert len(alleles) == 2 and alleles[0]["mutated_common_window"] != alleles[1]["mutated_common_window"]
    write_json("TYMSOS_11bp_deletion_equivalence_check.json", {"common_reference": dd["reference"],
        "start0": dd["start0"], "alleles": alleles, "same_resulting_local_sequence": False})
    write_json("evidence_audit.json", {"summary": summaries, "events": event_rows,
        "cross_GX_duplicate_UMIs_in_support_record_regional_representatives": 0,
        "all_source_cell_state_sums_reconciled": True, "no_new_coverage_threshold_defined": True})
    return events, summaries


def box(ax, x, y, w, h, text, size=13):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="#F6F8F9", edgecolor="#9BA6AD", lw=.85))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=size, linespacing=1.28)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color="#596A75", lw=1))


def mixscale_figure():
    fig, ax = plt.subplots(figsize=(16, 9))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    fig.suptitle("3 · Independent DEG-testing and scoring routes", y=.97, fontsize=21)
    box(ax, .12, .84, .76, .085,
        "Saved RNA counts and frozen sgRNA assignments\nPilot groups: TYMSOS 681 · UBA1 538 · NTC 19,342 records", 14)
    arrow(ax, (.30, .84), (.25, .785)); arrow(ax, (.7, .84), (.73, .785))
    ax.text(.25, .765, "SCEPTRE · per sgRNA", ha="center", fontsize=17, color=BLUE)
    ax.text(.74, .765, "Mixscale · per target", ha="center", fontsize=17, color=BLUE)
    ax.plot([.49, .49], [.10, .78], color="#DDDDDD", lw=.8, ls="--")
    box(ax, .06, .59, .38, .12, "One sgRNA group vs NTC\nRaw RNA counts", 15)
    arrow(ax, (.25, .59), (.25, .53))
    box(ax, .06, .40, .38, .13, "SCEPTRE testing\nPermutations · BH FDR < 0.1", 15)
    arrow(ax, (.25, .40), (.25, .34))
    box(ax, .06, .20, .38, .14, "Per-sgRNA DEG lists\nTYMSOS: 0, 5, 2, 2, 2\nUBA1: 1, 0, 1, 0, 2", 14)
    ax.text(.25, .135, "These DEG lists are not passed to Mixscale.", ha="center", fontsize=12)
    right = [
        (.645, .09, "Normalize RNA → 2,000 HVGs → shared PCA\nPCA fit: 81,201 records; use PCs 1–10", 13),
        (.485, .13, "Pool sgRNA-assigned records by target\n5 designed sgRNAs per target in this pilot\nSeurat DEG vs NTC: TYMSOS 1 (LINC02315), UBA1 1 (ZFP42)", 12.3),
        (.345, .105, "20 nearest NTC records per record (PCs 1–10)\nSubtract their mean expression for selected genes", 13),
        (.215, .095, "Scoring entry: original min.de.genes = 5 → not scored\nLowered to 1 → score computed", 13),
        (.095, .085, "Use the background-corrected DEG direction\nExploratory single-gene score", 14),
    ]
    for i, (y, h, text, size) in enumerate(right):
        box(ax, .525, y, .43, h, text, size)
        if i:
            arrow(ax, (.74, right[i-1][0]), (.74, y+h))
    fig.text(.53, .067, "DE screen: Wilcoxon; Bonferroni-adjusted P < 0.05.", fontsize=11)
    fig.text(.06, .025, "The minimum DEG count is a scoring setting, not an experiment-quality standard. No shared DEG list between routes.", fontsize=12)
    save(fig, "Fig3_Mixscale_scoring_limits", "两条独立路线。Mixscale先建立共享PCA，随后按靶点筛选DEG；筛选结果参与背景差值计算所用基因集合的构建。再通过前10个PC寻找20个最近邻NTC，计算表达差值，最后进入评分。两个靶点各仅一个DEG，降低最低基因数至1后得到探索性单基因分数。五条sgRNA是本次两个靶点的设计情况，并非Mixscale要求。", [REPO / "scripts/run_two_target_mixscale.R", REPO / "scripts/run_mixscale_lower_minimum.R", CATALOGUE, WORKBOOK])


def sequence_figure(guide, data):
    row = next(r for r in data["guides"] if r["guide"] == guide)
    display = next(r for r in data["display"] if r["guide"] == guide)
    variants = sorted([v for v in data["variants"] if v["guide"] == guide], key=lambda v: (v["pos0"], v["kind"], v["alt"]))
    reference, start = display["reference"], display["start0"]
    n = len(variants)
    fig, ax = plt.subplots(figsize=(16, max(4.1, 2.7+.49*n)))
    fig.subplots_adjust(left=.025, right=.99, top=.79, bottom=.24)
    ax.set(xlim=(-12, len(reference)+37), ylim=(-n-.4, 2))
    ax.axvspan(row["start0"]-start-.5, row["end0"]-start-.5, color="#edf0f1", zorder=0)
    for i, base in enumerate(reference):
        ax.text(i, 1, base, ha="center", va="center", family="monospace", fontsize=10, color="#555555")
    ax.text(-1.8, 1, "Reference", ha="right", va="center", fontsize=11)
    if row["cut0"] is not None:
        ax.axvline(row["cut0"]-start-.5, color=RED, lw=.8, ls="--")
        ax.text(row["cut0"]-start, 1.7, "Predicted cut", color=RED, ha="center", fontsize=10)
    ax.text(len(reference)+2, 1, "Target-guide RNA\nALT / informative\nUMI families", va="center", fontsize=10.3)
    ax.text(len(reference)+20, 1, "NTC RNA\nALT / informative\nUMI families", va="center", fontsize=10.3)
    for j, v in enumerate(variants):
        y, pos = -j, v["pos0"]-start
        ax.text(-1.8, y, f"{v['kind']} {v['length_change']:+d} bp" if v["kind"] != "SNV" else f"{v['ref']}>{v['alt']}", ha="right", va="center", fontsize=11)
        for i, base in enumerate(reference):
            ax.text(i, y, base, ha="center", va="center", family="monospace", fontsize=10, color="#b7b7b7")
        if v["kind"] == "DEL":
            a, b = pos+1, pos+len(v["ref"])
            ax.plot([max(0, a)-.4, min(len(reference), b)-.6], [y, y], color=RED, lw=2.2)
        elif v["kind"] == "INS":
            ax.text(pos+.5, y, "▼", color=RED, ha="center", va="center", fontsize=11)
            ax.text(pos+.5, y+.30, "+"+v["alt"][len(v["ref"]):], color=RED, ha="center", fontsize=10, family="monospace")
        else:
            ax.text(pos, y, v["alt"], color=RED, ha="center", va="center", fontsize=10, family="monospace")
        for prefix, x in (("guide", len(reference)+2), ("NTC", len(reference)+20)):
            total = sum(v[prefix+"_"+k] for k in ("alt", "ref", "other"))
            ax.text(x, y, f"{v[prefix+'_alt']} / {total}", va="center", fontsize=12)
    ax.axis("off")
    fig.suptitle(f"{guide}\n{row['chrom']}:{row['start0']+1:,}–{row['end0']:,}; sgRNA region shaded", fontsize=16, y=.975)
    fig.text(.035, .15, "UMI families: sample × cell barcode × corrected UMI × gene tag (GX). Counts are not reads or cell-record counts.", fontsize=10)
    fig.text(.035, .10, "Informative = ALT + REF + OTHER; conflicting families excluded. Each row is one normalized sgRNA–variant pair.", fontsize=10)
    cap = " Regional scan reached its read limit; counts describe the saved observations." if row["scan_capped"] else ""
    fig.text(.035, .05, "RNA sequence evidence only; not a DNA genotype or functional KO call." + cap, fontsize=10)
    save(fig, f"RNA_indel_sequence_{guide}", f"{guide}的RNA序列示意。ALT为支持该规范化事件的UMI-family数；informative为ALT、REF及OTHER三类之和，不含家族内部冲突。各行独立计数，同一sgRNA可出现多个不同事件。同一条记录可贡献多个不同UMI。", [SEQUENCE, REPO / "scripts/plot_editing_rna_cells.py"])


def pca_figures(labels):
    old = module("revision_shared_pca", OPS / "20261110_shared_PCA_preview/plot_shared_pca_indels.py")
    track(OPS / "20261109_target_indel_pca/target_coordinates.zip")
    track(PROJECT / "reports/20261029_two_figure_briefing/rna_indel_pca_20261109/render_receipt.json")
    atlas, cases, variants, _ = old.load_existing_labels()
    proof = read_json(OPS / "20261110_shared_PCA_preview/coordinate_receipt.json")
    payload = track(COORDS).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == proof["coordinates_sha256"]
    with np.load(io.BytesIO(payload), allow_pickle=False) as z:
        xy, guides, barcodes = [z[k] for k in ("xy", "guide", "cell_barcode")]
        ntc = xy[z["role"] == "non_targeting_control"]
        variance = 100*z["variance_ratio"]
    lo, hi = xy.min(0), xy.max(0)
    limits = np.array([lo-.045*(hi-lo), hi+.045*(hi-lo)]).T
    # Reuse the unchanged point-drawing function; intercept its save/close only
    # to add a definition footer and send PNG+SVG to the new directory.
    for name in ("TYMSOS", "UBA1"):
        chosen = sorted([v for v in variants if v["target_gene_id"] == cases[name]["id"] and v["kind"] in ("INS", "DEL") and v["guide_alt"] > 0], key=identity)
        mask = np.char.startswith(guides, name+"-")
        original_save, original_close = plt.Figure.savefig, plt.close
        held = []
        try:
            def capture(fig, *args, **kwargs):
                held.append(fig)
                original_save(fig, *args, **kwargs)
            plt.Figure.savefig = capture
            plt.close = lambda *a, **kw: None
            old.OUT = HERE
            old.plot_target(cases[name], chosen, xy[mask], guides[mask], barcodes[mask], ntc, limits, variance)
        finally:
            plt.Figure.savefig, plt.close = original_save, original_close
        assert len(held) == 1
        fig = held[0]
        fig.set_size_inches(12.75, 9.7)
        fig.subplots_adjust(bottom=.145, top=.855, hspace=.32)
        for ax in fig.axes:
            if ax.get_visible():
                ax.set_xticks([t for t in ax.get_xticks() if limits[0, 0] <= t <= limits[0, 1]])
                ax.set_yticks([t for t in ax.get_yticks() if limits[1, 0] <= t <= limits[1, 1]])
        fig._supxlabel.set_y(.083)
        fig._supylabel.set_x(.014)
        for legend in fig.legends:
            for text in legend.get_texts():
                if text.get_text() == "Supports shown indel":
                    text.set_text("RNA supports shown indel")
        fig.text(.07, .045, "Event n = unique (sgRNA, barcode) records with ≥1 ALT UMI family; n is not the number of RNA families.", fontsize=9)
        fig.text(.07, .019, "One sgRNA may support multiple distinct normalized indel events. RNA evidence does not establish DNA or functional KO.", fontsize=9)
        save(fig, f"{name}_indels_shared_PCA", f"{name}事件级PCA。完全复用已保存的共享PCA坐标和原来的事件支持记录；n逐事件按(sgRNA,barcode)去重，不是UMI-family数。TYMSOS各事件n=1、1、5、1、1，UBA1为8、3、7、13；红点仅表示对应RNA indel证据。", [LABELS, COORDS, OPS / "20261110_shared_PCA_preview/plot_shared_pca_indels.py"])


def summary_figures():
    old = module("revision_meeting_summary", REPO / "scripts/plot_editing_meeting_summary.py")
    support, effects = old.sequence_counts(SEQUENCE, CODING)
    summary = old.guide_sequence_summary(SEQUENCE, CODING)
    # Match the original embedded Nov 7 images to saved render receipts.
    with ZipFile(track(PPT)) as z:
        embedded = {n: hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith("ppt/media/")}
    matches = []
    for label, folder in (("events", "scientific_style_final"), ("guides", "guide197_final")):
        rec = read_json(OPS / f"20261107_meeting_figures/{folder}/render_receipt.json")
        h = rec["figures_sha256"]["sequence_changes.png"]
        assert h in embedded.values()
        matches.append({"unit": label, "embedded_media": next(k for k, v in embedded.items() if v == h), "source_receipt": str(OPS / f"20261107_meeting_figures/{folder}/render_receipt.json")})
    write_json("Nov7_PPT_figure_matches.json", matches)
    fig = old.sequence_figure(support, effects, percent=False)
    fig.set_size_inches(12, 6.6)
    fig.subplots_adjust(bottom=.29, top=.89)
    fig.axes[0].set_title("RNA variants (counted separately for each sgRNA)")
    fig.axes[0].set_ylabel("Number of sgRNA–variant pairs")
    fig.axes[1].set_title("Predicted indel effects")
    fig.text(.09, .08, "One pair = one sgRNA + one normalized RNA variant. Same variant under different sgRNAs counts separately.", fontsize=9.5)
    fig.text(.09, .035, "Family-support categories count UMI families, not reads or cells. Predicted RNA effects do not establish functional KO.", fontsize=9.5)
    save(fig, "Nov7_RNA_indel_event_summary", "11月7日原事件汇总图的单位澄清版。保留替换100、插入43、缺失102项及所有原柱高；一项为一条sgRNA与一个规范化变异的组合，不是细胞数。同一变异出现在不同sgRNA组时分别计数。颜色表示支持该组合的UMI-family数。", [SEQUENCE, CODING, PPT, REPO / "scripts/plot_editing_meeting_summary.py"])
    fig = old.guide_sequence_figure(summary, percent=True)
    fig.set_size_inches(12, 6.3)
    fig.subplots_adjust(bottom=.25, top=.87)
    fig.text(.09, .09, "Denominator: 197 selected sgRNAs; each sgRNA is counted once within a category. Categories can overlap.", fontsize=10)
    fig.text(.09, .04, "A single sgRNA can support multiple indel events. RNA sequence evidence is not a DNA-editing or functional-KO rate.", fontsize=10)
    save(fig, "Nov7_RNA_indel_sgRNA_summary", "11月7日原sgRNA比例图的图注补充版。分母为197条已选sgRNA，每类内按sgRNA去重；114条有任一RNA变化，101条有indel，81条有预计移码RNA，各类别可以重叠。同一sgRNA可以支持多个indel事件。这不是DNA编辑率。", [SEQUENCE, CODING, PPT, REPO / "scripts/plot_editing_meeting_summary.py"])


def documents(design, coverage):
    captions = "\n\n".join(f"{r['stem']}\n{r['caption']}" for r in FIGURES)
    (HERE / "中文图注.txt").write_text(captions+"\n", encoding="utf-8")
    table("绘图索引.csv", [{"figure": f["stem"], "png": str(HERE / (f["stem"]+".png")),
        "editable_svg": str(HERE / (f["stem"]+".svg")), "source_files": ";".join(f["sources"]),
        "redraw_script": str(Path(__file__).resolve())} for f in FIGURES])
    html = "<!doctype html><meta charset='utf-8'><title>RNA证据核查与修图</title><style>body{font:17px 'Microsoft Yahei',sans-serif;max-width:1250px;margin:30px auto;line-height:1.6}img{max-width:100%;border:1px solid #ddd}section{margin:40px 0}a{color:#286281}</style><h1>RNA证据核查与修图</h1><p><a href='中文核查结果.md'>中文核查结果</a> · <a href='indel_supported_record_evidence_40.csv'>40条支持记录明细</a> · <a href='event_counts_and_PCA_n.csv'>事件与n核对</a></p>"
    for f in FIGURES:
        html += f"<section><h2>{f['stem']}</h2><a href='{f['stem']}.svg'>可编辑SVG</a><img src='{f['stem']}.png'><p>{f['caption']}</p></section>"
    (HERE / "index.html").write_text(html, encoding="utf-8")
    thumbs = []
    for f in FIGURES:
        with Image.open(HERE / (f["stem"]+".png")) as im:
            thumb = ImageOps.contain(im.convert("RGB"), (800, 500))
        tile = Image.new("RGB", (830, 550), "white")
        tile.paste(thumb, ((830-thumb.width)//2, 35))
        ImageDraw.Draw(tile).text((15, 12), f["stem"], fill="black")
        thumbs.append(tile)
    contact = Image.new("RGB", (1660, 550*((len(thumbs)+1)//2)), "#eeeeee")
    for i, im in enumerate(thumbs):
        contact.paste(im, ((i%2)*830, (i//2)*550))
    contact.save(HERE / "all_revisions_preview.png")


def main():
    assert HERE.name == "20261114_rna_evidence_revision"
    sys.path.insert(0, str(REPO / "scripts"))
    style()
    version_checks = []
    for name, functions in (("analyze_guided_rna_alleles.py", ("Read", "raw_events", "covers_reference", "normalize_events")),
                            ("explore_editing_rna_cells.py", ("family_state", "load_cache", "extract_region"))):
        past = subprocess.check_output(["git", "show", "09b6ca35ccce371988f6e14a49a1f16732257ae9:scripts/"+name], cwd=REPO, text=True, encoding="utf-8")
        current = track(REPO / "scripts" / name).read_text(encoding="utf-8-sig")
        def nodes(text):
            return {x.name: ast.dump(x, include_attributes=False) for x in ast.parse(text).body if isinstance(x, (ast.ClassDef, ast.FunctionDef))}
        left, right = nodes(past), nodes(current)
        for function in functions:
            assert left[function] == right[function], (name, function)
            version_checks.append({"script": name, "function": function, "matches_recorded_analysis_commit": True})
    write_json("historical_code_verification.json", {"recorded_analysis_commit": "09b6ca35ccce371988f6e14a49a1f16732257ae9", "checks": version_checks})
    # All old figure exports and the formal PPT are protected by hash checks.
    for path in (REPO / "meeting_figures/20261114_perturbseq").glob("Fig*.*"):
        if path.suffix in (".png", ".svg"):
            track(path)
    for p in ("analyze_guided_rna_alleles.py", "explore_editing_rna_cells.py", "cache_guided_rna_alignments.py",
              "plot_editing_rna_cells.py", "plot_indel_pca_by_target.py", "run_two_target_mixscale.R", "run_mixscale_lower_minimum.R"):
        track(REPO / "scripts" / p)
    data = read_json(SEQUENCE)["company"]
    labels = read_json(LABELS)["company"]
    design = design_audit()
    events, coverage = evidence_audit(data, labels)
    assert len(events) == 9
    mixscale_figure()
    for guide in ("TYMSOS-2", "TYMSOS-4", "TYMSOS-5", "UBA1-2", "UBA1-3"):
        sequence_figure(guide, data)
    pca_figures(labels)
    summary_figures()
    documents(design, coverage)
    assert len(FIGURES) == 10
    assert all(sha(path) == v["sha256"] for path, v in INPUTS.items()), "Source modified during work"
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    write_json("source_manifest.json", {"code_commit_at_run": commit, "script_sha256": sha(__file__),
        "inputs": INPUTS, "all_tracked_inputs_unchanged": True,
        "no_BAM_or_expression_matrix_read": True, "no_analysis_refit": True})
    write_json("render_checks.json", FIGURES)
    print("AUDIT_OK", coverage, flush=True)


if __name__ == "__main__":
    main()
