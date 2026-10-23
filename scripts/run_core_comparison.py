"""Compare accepted per-chip core results without rereading counts or FASTQs."""
import argparse
from datetime import datetime
import itertools
import json
import os
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

from manage_core_study import filesystem_kind
from perturbseq_core_io import atomic_json, identity, revision, sha256


ROLES = ("non_targeting_control", "intergenic_control", "positive_control", "targeting")
PAYLOADS = ("summary.json", "tables/control_annotations.csv",
            "tables/guide_effect_metadata.csv", "tables/rna_effect_profiles.npz",
            "tables/target_support.csv")


def read_json(path):
    if path.stat().st_size > 8 * 1024 * 1024:
        raise ValueError("JSON metadata budget exceeded")
    return json.loads(path.read_text(encoding="utf8"))


def read_bundle(spec, smoke=False):
    root = Path(spec["path"]).resolve()
    if not Path(spec["path"]).is_absolute():
        raise ValueError("Absolute input paths required")
    manifest_path = root / "output_manifest.json"
    if not re.fullmatch(r"[0-9a-f]{64}", spec["manifest_sha256"]):
        raise ValueError("Accepted input manifest SHA256 required")
    if sha256(manifest_path) != spec["manifest_sha256"]:
        raise ValueError("Input manifest differs")
    if (root / "owner.lock").exists() or read_json(root / "status.json").get("stage") != "complete":
        raise ValueError("Input is not complete or still has an owner")
    manifest = read_json(manifest_path)
    snapshots = {}
    for rel in PAYLOADS:
        path = (root / rel).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Input payload escapes its run")
        proof = manifest[rel]
        snapshots[rel] = identity(path)
        if path.stat().st_size != proof["bytes"] or sha256(path) != proof["sha256"]:
            raise ValueError("Input payload differs: " + rel)
    summary = read_json(root / "summary.json")
    if summary.get("input_scope", "independent_capture") != "independent_capture":
        raise ValueError("Merged diagnostic records cannot be used as an independent chip")
    if summary.get("scientific_result") is not True and not smoke:
        raise ValueError("Engineering input cannot produce a scientific comparison")
    if summary.get("chip_id", spec["chip_id"]) != spec["chip_id"]:
        raise ValueError("Chip label differs from source")
    counts_sha = summary["source_identity"]["sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", counts_sha):
        raise ValueError("Raw-count source identity missing")
    legacy = "core_reference" in summary
    if legacy and summary.get("code_commit") != "c82bcf418731122de2f3244a4bbb1baed6225d7d":
        raise ValueError("Unrecognized historical summary; explicit adapter needed")
    scale = summary["core_reference"]["normalization_target_umi"] if legacy else summary["normalization_target_umi"]
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid normalization scale")
    matched_cells = summary["matched_cells"]
    if isinstance(matched_cells, bool) or int(matched_cells) != matched_cells or matched_cells < 1:
        raise ValueError("Invalid matched cell count")
    ntc_cells = int(summary["core_reference"]["all_pure_ntc_cells"] if legacy else summary["reference_cells"])
    catalogue = pd.read_csv(root / "tables/control_annotations.csv", keep_default_na=False)
    required = ["guide", "target_gene_name", "source_control_type"]
    catalogue = catalogue[required].copy()
    if catalogue.guide.duplicated().any() or not set(catalogue.source_control_type).issubset(ROLES):
        raise ValueError("Duplicate guide or unknown control role")
    metadata = pd.read_csv(root / "tables/guide_effect_metadata.csv", keep_default_na=False)
    if metadata.guide.duplicated().any() or not metadata.reference_overlap_cells.eq(0).all():
        raise ValueError("Duplicate guide or overlapping NTC reference")
    profile_path = root / "tables/rna_effect_profiles.npz"
    with zipfile.ZipFile(profile_path) as archive:
        if sum(i.file_size for i in archive.infolist()) > 2 * 1024**3:
            raise ValueError("Profile expansion exceeds comparison memory budget")
    with np.load(profile_path, allow_pickle=False) as data:
        guides, genes = data["guides"].astype(str), data["gene_ids"].astype(str)
        full, matched = data["log2fc"].copy(), data["matched_first"].copy()
        pseudocount = float(data["pseudocount"])
    if len(set(guides)) != len(guides) or len(set(genes)) != len(genes):
        raise ValueError("Duplicate profile identities")
    if set(guides) != set(metadata.guide) or not set(guides).issubset(set(catalogue.guide)):
        raise ValueError("Profile/guide table identity mismatch")
    if full.shape != (len(guides), len(genes)) or matched.shape != full.shape:
        raise ValueError("Profile dimensions differ")
    metadata = metadata.set_index("guide").loc[guides].reset_index()
    if "direct_feature_index" not in metadata or pseudocount != 1:
        raise ValueError("Direct-target exclusion provenance or pseudocount missing")
    excluded = np.zeros(full.shape, dtype=bool)
    for row, value in enumerate(metadata.direct_feature_index):
        if value == "" or pd.isna(value):
            continue
        index = float(value)
        if not np.isfinite(index) or not index.is_integer() or not 0 <= index < len(genes):
            raise ValueError("Invalid direct-target feature index")
        excluded[row, int(index)] = True
    if (np.isinf(full).any() or np.isinf(matched).any() or
        not np.array_equal(np.isnan(full), excluded) or
        not np.array_equal(np.isnan(matched), excluded)):
        raise ValueError("Profile missing values differ from declared direct-target exclusions")
    labels = catalogue.set_index("guide").loc[guides]
    if (metadata.source_control_type.to_numpy() != labels.source_control_type.to_numpy()).any():
        raise ValueError("Profile control role differs from catalogue")
    if (metadata.target.to_numpy() != labels.target_gene_name.to_numpy()).any():
        raise ValueError("Profile target differs from catalogue")
    expected_reference = ntc_cells - np.where(
        metadata.source_control_type.eq("non_targeting_control"), metadata.cells, 0)
    if (expected_reference < 2).any() or not np.array_equal(
            metadata.reference_cells.to_numpy(), expected_reference):
        raise ValueError("Reference counts do not implement whole-guide NTC leave-out")
    if (metadata.cells < summary["matched_cells"]).any():
        raise ValueError("A profile lacks the declared matched cell count")
    source = dict(chip_id=spec["chip_id"], root=root, manifest_sha256=spec["manifest_sha256"],
                  catalogue=catalogue, metadata=metadata, genes=genes, guides=guides,
                  full=full, matched=matched, scale=float(scale),
                  matched_cells=int(summary["matched_cells"]), snapshots=snapshots,
                  source_cells=summary["source_identity"]["source_shape"][0] if legacy else summary["source_cells"],
                  analysis_cells=summary["source_identity"]["selected_cells"] if legacy else summary["analysis_cells"],
                  ntc_cells=ntc_cells,
                  source_counts_sha256=counts_sha,
                  source_commit=summary["code_commit"] if legacy else summary["code"]["commit"],
                  target_support=pd.read_csv(root / "tables/target_support.csv"))
    return source


def compare(bundles, minimum_common_features=2):
    """Align by actual RNA/guide IDs; never average unmatched gene columns."""
    if len(bundles) < 2 or len({b["chip_id"] for b in bundles}) != len(bundles):
        raise ValueError("At least two distinct chips required")
    if len({b["source_counts_sha256"] for b in bundles}) != len(bundles):
        raise ValueError("The same count source cannot represent multiple chips")
    if len({b["scale"] for b in bundles}) != 1 or len({b["matched_cells"] for b in bundles}) != 1:
        raise ValueError("Normalization depth and matched cell count must agree")
    canonical = bundles[0]["catalogue"].sort_values("guide").reset_index(drop=True)
    for b in bundles[1:]:
        if not canonical.equals(b["catalogue"].sort_values("guide").reset_index(drop=True)):
            raise ValueError("Guide targets/control roles differ between chips")
    common = sorted(set.intersection(*(set(b["genes"]) for b in bundles)))
    if len(common) < max(2, minimum_common_features):
        raise ValueError("Insufficient common RNA feature space")
    guide_metrics, support = [], []
    projected = {}
    for b in bundles:
        lookup = {g: i for i, g in enumerate(b["genes"])}
        columns = [lookup[g] for g in common]
        projected[b["chip_id"]] = b["full"][:, columns]
        table = b["metadata"].copy()
        table.insert(0, "chip_id", b["chip_id"])
        table["common_features"] = len(common)
        # One already-frozen matched draw per guide, identical n and RNA space.
        table["effective_features"] = np.isfinite(b["matched"][:, columns]).sum(axis=1)
        table["matched_rms_log2fc"] = np.sqrt(np.nanmean(b["matched"][:, columns] ** 2, axis=1))
        guide_metrics.append(table)
        support.append(b["target_support"].assign(chip_id=b["chip_id"]))
    pairs = []
    for left, right in itertools.combinations(bundles, 2):
        li = {g: i for i, g in enumerate(left["guides"])}
        ri = {g: i for i, g in enumerate(right["guides"])}
        labels = canonical.set_index("guide")
        for guide in sorted(li.keys() & ri.keys()):
            a, b = projected[left["chip_id"]][li[guide]], projected[right["chip_id"]][ri[guide]]
            keep = np.isfinite(a) & np.isfinite(b)
            a, b = a[keep], b[keep]
            value = float(np.corrcoef(a, b)[0, 1]) if len(a) > 1 and np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan
            pairs.append(dict(chip_a=left["chip_id"], chip_b=right["chip_id"], guide=guide,
                              target=labels.loc[guide, "target_gene_name"],
                              source_control_type=labels.loc[guide, "source_control_type"],
                              common_features=len(common), compared_features=len(a), pearson=value))
    columns = ["chip_a", "chip_b", "guide", "target", "source_control_type", "common_features", "compared_features", "pearson"]
    return common, pd.concat(guide_metrics, ignore_index=True), pd.DataFrame(pairs, columns=columns), pd.concat(support, ignore_index=True)


def figures(out, bundles, metrics, pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chips = [b["chip_id"] for b in bundles]
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    x = np.arange(len(chips))
    ax.bar(x-.18, [b["source_cells"] for b in bundles], .36, label="Input cells")
    ax.bar(x+.18, [b["analysis_cells"] for b in bundles], .36, label="RNA-QC single-guide cells")
    ax.set(xticks=x, xticklabels=chips, ylabel="Cells", title="Cell coverage by capture")
    ax.legend()
    fig.savefig(out/"figures/01_cell_coverage.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4), layout="constrained")
    for role, offset, color in [("non_targeting_control", -.17, "#337ca0"), ("targeting", .17, "#bf6d38")]:
        for i, chip in enumerate(chips):
            values = metrics.loc[metrics.chip_id.eq(chip) & metrics.source_control_type.eq(role), "matched_rms_log2fc"]
            if len(values):
                box = ax.boxplot([values], positions=[i+offset], widths=.28, patch_artist=True,
                                 manage_ticks=False, medianprops={"color":color},
                                 whiskerprops={"color":color}, capprops={"color":color})
                box["boxes"][0].set_facecolor(color)
        ax.plot([], [], color=color, linewidth=6, label="NTC" if role.startswith("non_") else "Targeting")
    ax.set(xticks=x, xticklabels=chips, ylabel="RMS RNA log2FC", title="Same RNA features and cell count; one frozen matched draw")
    ax.legend()
    fig.savefig(out/"figures/02_matched_background.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for role, ax, title in zip(["non_targeting_control", "targeting"], axes, ["NTC", "Targeting"]):
        values = np.full((len(chips), len(chips)), np.nan)
        for i, a in enumerate(chips):
            for j, b in enumerate(chips):
                if i == j:
                    continue
                subset = pairs.loc[pairs.source_control_type.eq(role) &
                    ((pairs.chip_a.eq(a) & pairs.chip_b.eq(b)) | (pairs.chip_a.eq(b) & pairs.chip_b.eq(a))), "pearson"].dropna()
                if len(subset):
                    values[i, j] = subset.median()
                    ax.text(j, i, f"{values[i,j]:.2f}\nn={len(subset)}", ha="center", va="center", fontsize=8)
        cmap = plt.get_cmap("coolwarm").with_extremes(bad="#e5e7eb")
        im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=-1, vmax=1)
        ax.set(xticks=x, yticks=x, xticklabels=chips, yticklabels=chips, title=title+" same-guide agreement")
        ax.tick_params(axis="x", rotation=45)
    fig.colorbar(im, ax=axes, label="Median Pearson r across common eligible guides")
    fig.savefig(out/"figures/03_cross_chip_agreement.png", dpi=160)
    plt.close(fig)


def verify(out):
    out = Path(out).resolve()
    if (out/"owner.lock").exists() or read_json(out/"status.json").get("stage") != "complete":
        raise ValueError("Comparison incomplete or owned")
    manifest = read_json(out/"output_manifest.json")
    actual = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()} - {"status.json", "output_manifest.json"}
    if actual != set(manifest):
        raise ValueError("Output coverage differs")
    for rel, proof in manifest.items():
        p = (out/rel).resolve()
        if not p.is_relative_to(out) or p.stat().st_size != proof["bytes"] or sha256(p) != proof["sha256"]:
            raise ValueError("Output differs: "+rel)
    for rel in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", (out/"总体与芯片间对照报告.md").read_text(encoding="utf8")):
        if rel not in manifest:
            raise ValueError("Report image absent from manifest")
    return dict(complete=True, files=len(manifest), scientific_result=read_json(out/"summary.json")["scientific_result"])


def run(config_path, smoke=False):
    config = read_json(Path(config_path))
    if config.get("schema_version") != 1 or config.get("capture_units_confirmed") is not True:
        raise ValueError("Confirmed capture design and comparison schema required")
    expected = config["expected_chips"]
    if len(expected) < 2 or len(set(expected)) != len(expected) or set(expected) != {s["chip_id"] for s in config["inputs"]} or len(expected) != len(config["inputs"]):
        raise ValueError("Exact requested chip coverage required")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,95}", c) for c in expected):
        raise ValueError("Unsafe chip label")
    code = revision(config["expected_commit"], smoke)
    out = Path(config["output"]).resolve()
    parent = Path(config["allowed_output_parent"]).resolve()
    repo = Path(__file__).resolve().parents[1]
    roots = [Path(s["path"]).resolve() for s in config["inputs"]]
    if len(set(roots)) != len(roots):
        raise ValueError("Duplicate input run directories")
    if not Path(config["output"]).is_absolute() or out.parent != parent or out.is_relative_to(repo) or any(out.is_relative_to(p) or p.is_relative_to(out) for p in roots):
        raise ValueError("Separate external output required")
    if not smoke:
        spec = config["comparability_review"]
        p = Path(spec["path"]).resolve()
        if sha256(p) != spec["sha256"]:
            raise ValueError("Comparison review differs")
        review = read_json(p)
        if review.get("compatible") is not True or not str(review.get("basis", "")).strip() or review.get("input_manifests") != {s["chip_id"]:s["manifest_sha256"] for s in config["inputs"]}:
            raise ValueError("Count/reference compatibility review is missing or stale")
    out.mkdir(exist_ok=False)
    owner = dict(pid=os.getpid(), started_at=datetime.now().astimezone().isoformat(), code=code)
    atomic_json(out/"owner.lock", owner)
    try:
        atomic_json(out/"status.json", dict(stage="reading_inputs", **owner))
        atomic_json(out/"config.resolved.json", config)
        if not smoke:
            atomic_json(out/"comparability_review.json", review)
        specs = {s["chip_id"]:s for s in config["inputs"]}
        bundles = [read_bundle(specs[chip], smoke) for chip in expected]
        common, metrics, pairs, support = compare(bundles, config["minimum_common_features"])
        (out/"tables").mkdir()
        (out/"figures").mkdir()
        pd.DataFrame({"gene_id":common}).to_csv(out/"tables/common_response_genes.csv", index=False)
        metrics.to_csv(out/"tables/chip_guide_metrics.csv", index=False)
        pairs.to_csv(out/"tables/cross_chip_guide_agreement.csv", index=False)
        support.to_csv(out/"tables/target_support_by_chip.csv", index=False)
        bundles[0]["catalogue"].to_csv(out/"tables/guide_catalogue.csv", index=False)
        columns = ["chip_id", "source_cells", "analysis_cells", "ntc_cells", "source_commit",
                   "source_counts_sha256", "manifest_sha256"]
        pd.DataFrame([{k:b[k] for k in columns} for b in bundles]).to_csv(out/"tables/chip_summary.csv", index=False)
        figures(out, bundles, metrics, pairs)
        audits = []
        for b in bundles:
            if (b["root"]/"owner.lock").exists() or read_json(b["root"]/"status.json").get("stage") != "complete":
                raise ValueError("Source ownership/state changed during comparison")
            if sha256(b["root"]/"output_manifest.json") != b["manifest_sha256"]:
                raise ValueError("Source manifest changed during comparison")
            for rel, before in b["snapshots"].items():
                after = identity(b["root"]/rel)
                if before != after:
                    if filesystem_kind(b["root"]) in {"cifs", "smb3"} and all(before[k] == after[k] for k in ["bytes", "mtime_ns", "device"]):
                        audits.append(dict(chip_id=b["chip_id"], file=rel, before=before, after=after))
                    else:
                        raise ValueError("Source changed during comparison: "+rel)
        if not smoke and sha256(Path(config["comparability_review"]["path"])) != config["comparability_review"]["sha256"]:
            raise ValueError("Comparability review changed during comparison")
        summary = dict(scientific_result=not smoke, chip_ids=expected, code=code,
            common_features=len(common), matched_cells=bundles[0]["matched_cells"],
            normalization_target_umi=bundles[0]["scale"], comparison_unit="capture and guide",
            raw_data_read=False, counts_read=False, cells_pooled=False, biological_replicates_established=False,
            validated_functional_hits=False, inferential_tests_run=False, target_inode_audits=audits)
        atomic_json(out/"summary.json", summary)
        report = f"""# 总体与芯片间对照报告

{'这是合成/工程检查结果，不能用于科学汇报。' if smoke else '本报告比较已验收的各芯片结果，检查覆盖情况、阴性对照背景和同一sgRNA的反应是否一致。'}

本次纳入 {len(bundles)} 张芯片：{', '.join(expected)}。原始计数和各芯片报告保留原样。

## 各芯片提供了多少可用细胞

![各芯片细胞覆盖](figures/01_cell_coverage.png)

左侧柱表示输入细胞数，右侧表示通过RNA质控且得到单一sgRNA归属的细胞数。数量差异用于检查覆盖情况；不能仅凭柱子的高低判断实验成功或失败。逐靶点的可评价/覆盖不足状态保留在 tables/target_support_by_chip.csv。

## 阴性对照和targeting的RNA变化幅度

![等细胞数对照](figures/02_matched_background.png)

比较使用各芯片共有的 {len(common):,} 个RNA条目，每条sgRNA采用原分析已固定的一次 {summary['matched_cells']} 细胞抽样。纵轴是这些RNA条目log2FC的均方根，表示变化幅度而非变化方向；每条sgRNA贡献一个数值。NTC参考在各芯片内部建立，并排除被比较NTC guide自身的全部细胞。箱线图接近不等于已经证明没有扰动效应，偏离也不等于敲除成功。

这是重新对齐到共同RNA范围的描述性比较，不能把它与各芯片原来使用不同特征范围的数值直接混用。沿用原分析对该guide直接靶点RNA的排除，实际参与计算的条目数逐guide保存。共有RNA清单见 tables/common_response_genes.csv；它不代表全部转录组。某芯片某组没有箱线表示没有达到覆盖要求的guide，不能解释为变化为零。

## 同一条sgRNA在不同芯片上是否表现一致

![跨芯片同guide一致性](figures/03_cross_chip_agreement.png)

对两张芯片都达到覆盖要求的同一条sgRNA，使用全部合格细胞的平均RNA log2FC计算相关性；每格显示相关系数的中位数和可计算的guide数量。正相关表示共同RNA范围内的变化模式较一致，接近零表示一致性弱；灰色格表示未进行或无法计算的比较，包括对角线上的同一芯片。完整逐guide结果见 tables/cross_chip_guide_agreement.csv。

这些芯片来自同一实验的独立捕获，不能把芯片数当作独立生物学重复数。报告没有新增显著性检验、功能hit或细胞池化分析；低覆盖和未评价的靶点仍保留在表中。计数方法和参考定义的可比性由输入验收及随配置保存的审查记录约束。
"""
        (out/"总体与芯片间对照报告.md").write_text(report, encoding="utf8")
        manifest = {p.relative_to(out).as_posix():dict(bytes=p.stat().st_size, sha256=sha256(p))
                    for p in sorted(out.rglob("*")) if p.is_file() and p.name not in {"owner.lock", "status.json", "output_manifest.json"}}
        atomic_json(out/"output_manifest.json", manifest)
        atomic_json(out/"status.json", dict(stage="complete", **owner))
    except BaseException as error:
        atomic_json(out/"status.json", dict(stage="failed", error=repr(error), **owner))
        raise
    finally:
        (out/"owner.lock").unlink()
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    execute = commands.add_parser("run")
    execute.add_argument("--config", type=Path, required=True)
    execute.add_argument("--smoke", action="store_true")
    check = commands.add_parser("verify")
    check.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.smoke) if args.command == "run" else verify(args.output), ensure_ascii=False))
