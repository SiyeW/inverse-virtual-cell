import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_core_comparison as comparison


def bundle(root, chip, order=(0, 1, 2), legacy=False, scientific=True):
    root.mkdir()
    (root / "tables").mkdir()
    genes = np.array(["g1", "g2", "g3"])[list(order)]
    guides = np.array(["n1", "n2", "a1", "a2"])
    roles = ["non_targeting_control"] * 2 + ["targeting"] * 2
    targets = ["NTC", "NTC", "A", "A"]
    full = np.array([[.1, -.1, 0], [0, .1, -.1], [1, 2, 3], [2, 4, 6]])[:, list(order)]
    pd.DataFrame(dict(guide=guides, target_gene_name=targets, source_control_type=roles)).to_csv(root / "tables/control_annotations.csv", index=False)
    pd.DataFrame(dict(guide=guides[::-1], target=targets[::-1], source_control_type=roles[::-1],
                      cells=[20]*4, reference_cells=[40,40,20,20], reference_overlap_cells=[0]*4,
                      direct_feature_index=[""]*4)).to_csv(root / "tables/guide_effect_metadata.csv", index=False)
    pd.DataFrame([dict(target="A", source_control_type="targeting", eligible_guides=2,
                       assessment="descriptive_assessable")]).to_csv(root / "tables/target_support.csv", index=False)
    np.savez_compressed(root / "tables/rna_effect_profiles.npz", guides=guides, gene_ids=genes,
                        log2fc=full, matched_first=full, pseudocount=1.)
    summary = dict(chip_id=chip, scientific_result=scientific, normalization_target_umi=1000.,
                   matched_cells=10, source_cells=100, analysis_cells=80, reference_cells=40,
                   code=dict(commit="a"*40), source_identity=dict(sha256=hashlib.sha256(chip.encode()).hexdigest()))
    if legacy:
        summary = dict(scientific_result=scientific, matched_cells=10,
            code_commit="c82bcf418731122de2f3244a4bbb1baed6225d7d",
            source_identity=dict(source_shape=[100, 3], selected_cells=80,
                                 sha256=hashlib.sha256(chip.encode()).hexdigest()),
            core_reference=dict(normalization_target_umi=1000., all_pure_ntc_cells=40))
    (root / "summary.json").write_text(json.dumps(summary))
    (root / "status.json").write_text('{"stage":"complete"}')
    return renew(root, chip)


def renew(root, chip):
    manifest = {rel:dict(bytes=(root/rel).stat().st_size, sha256=comparison.sha256(root/rel)) for rel in comparison.PAYLOADS}
    (root / "output_manifest.json").write_text(json.dumps(manifest))
    return dict(chip_id=chip, path=str(root), manifest_sha256=comparison.sha256(root / "output_manifest.json"))


def config(tmp_path):
    specs = [bundle(tmp_path / "source1", "chip01", legacy=True),
             bundle(tmp_path / "source2", "chip02", order=(2, 0, 1))]
    value = dict(schema_version=1, capture_units_confirmed=True, expected_chips=["chip01", "chip02"],
        expected_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(comparison.__file__).parents[1], text=True).strip(),
        output=str(tmp_path / "comparison"), allowed_output_parent=str(tmp_path),
        inputs=specs, minimum_common_features=2)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(value))
    return path, value


def test_align_gene_and_guide_order_and_legacy_adapter(tmp_path):
    _, cfg = config(tmp_path)
    sources = [comparison.read_bundle(s) for s in cfg["inputs"]]
    genes, metrics, pairs, support = comparison.compare(sources)
    assert genes == ["g1", "g2", "g3"]
    assert np.allclose(pairs.pearson, 1)
    assert metrics.groupby("chip_id").guide.count().to_dict() == {"chip01":4, "chip02":4}
    assert len(support) == 2


@pytest.mark.parametrize("key,value", [("scale", 2000.), ("matched_cells", 11)])
def test_incompatible_comparison_units_stop(tmp_path, key, value):
    _, cfg = config(tmp_path)
    sources = [comparison.read_bundle(s) for s in cfg["inputs"]]
    sources[1][key] = value
    with pytest.raises(ValueError, match="must agree"):
        comparison.compare(sources)


def test_control_role_conflict_and_insufficient_common_space_stop(tmp_path):
    _, cfg = config(tmp_path)
    sources = [comparison.read_bundle(s) for s in cfg["inputs"]]
    bad = copy.deepcopy(sources)
    bad[1]["catalogue"].loc[0, "source_control_type"] = "targeting"
    with pytest.raises(ValueError, match="control roles"):
        comparison.compare(bad)
    sources[1]["genes"] = np.array(["other1", "other2", "other3"])
    with pytest.raises(ValueError, match="common RNA"):
        comparison.compare(sources)


def test_constant_profiles_are_unassessed_not_zero_correlation(tmp_path):
    _, cfg = config(tmp_path)
    sources = [comparison.read_bundle(s) for s in cfg["inputs"]]
    sources[1]["full"][0] = 0
    _, _, pairs, _ = comparison.compare(sources)
    assert pairs.loc[pairs.guide.eq("n1"), "pearson"].isna().all()


def test_duplicate_count_source_is_not_another_chip(tmp_path):
    _, cfg = config(tmp_path)
    sources = [comparison.read_bundle(s) for s in cfg["inputs"]]
    sources[1]["source_counts_sha256"] = sources[0]["source_counts_sha256"]
    with pytest.raises(ValueError, match="same count source"):
        comparison.compare(sources)


def test_duplicate_run_directory_rejected(tmp_path):
    path, cfg = config(tmp_path)
    cfg["inputs"][1]["path"] = cfg["inputs"][0]["path"]
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="Duplicate input"):
        comparison.run(path, smoke=True)


def test_changed_input_and_active_owner_stop(tmp_path):
    spec = bundle(tmp_path / "source", "chip01")
    (Path(spec["path"]) / "owner.lock").write_text("owned fixture")
    with pytest.raises(ValueError, match="owner"):
        comparison.read_bundle(spec)
    (Path(spec["path"]) / "owner.lock").unlink()
    (Path(spec["path"]) / "summary.json").write_text("{}")
    with pytest.raises(ValueError, match="payload differs"):
        comparison.read_bundle(spec)


def test_reference_overlap_rejected_even_with_valid_manifest(tmp_path):
    spec = bundle(tmp_path / "source", "chip01")
    path = Path(spec["path"]) / "tables/guide_effect_metadata.csv"
    table = pd.read_csv(path)
    table.loc[0, "reference_overlap_cells"] = 1
    table.to_csv(path, index=False)
    spec = renew(Path(spec["path"]), "chip01")
    with pytest.raises(ValueError, match="overlapping"):
        comparison.read_bundle(spec)


def test_wrong_reference_population_is_rejected(tmp_path):
    spec = bundle(tmp_path / "source", "chip01")
    root = Path(spec["path"])
    path = root / "tables/guide_effect_metadata.csv"
    table = pd.read_csv(path)
    table.loc[table.guide.eq("n1"), "reference_cells"] = 40
    table.to_csv(path, index=False)
    with pytest.raises(ValueError, match="whole-guide NTC"):
        comparison.read_bundle(renew(root, "chip01"))


def test_engineering_source_cannot_be_promoted(tmp_path):
    spec = bundle(tmp_path / "source", "chip01", scientific=False)
    with pytest.raises(ValueError, match="Engineering input"):
        comparison.read_bundle(spec)
    assert comparison.read_bundle(spec, smoke=True)["chip_id"] == "chip01"


def test_declared_direct_target_mask_is_preserved(tmp_path):
    _, cfg = config(tmp_path)
    spec = cfg["inputs"][0]
    root = Path(spec["path"])
    path = root / "tables/guide_effect_metadata.csv"
    table = pd.read_csv(path, keep_default_na=False)
    table.loc[table.guide.eq("a1"), "direct_feature_index"] = "0"
    table.to_csv(path, index=False)
    path = root / "tables/rna_effect_profiles.npz"
    with np.load(path) as z:
        arrays = {k:z[k].copy() for k in z.files}
    row = list(arrays["guides"]).index("a1")
    arrays["log2fc"][row, 0] = np.nan
    arrays["matched_first"][row, 0] = np.nan
    np.savez_compressed(path, **arrays)
    left = comparison.read_bundle(renew(root, "chip01"))
    right = comparison.read_bundle(cfg["inputs"][1])
    _, metrics, pairs, _ = comparison.compare([left, right])
    assert pairs.loc[pairs.guide.eq("a1"), "compared_features"].item() == 2
    assert np.isclose(pairs.loc[pairs.guide.eq("a1"), "pearson"].item(), 1)
    assert metrics.loc[metrics.chip_id.eq("chip01") & metrics.guide.eq("a1"), "effective_features"].item() == 2


def test_unexplained_nan_rejected(tmp_path):
    spec = bundle(tmp_path / "source", "chip01")
    root = Path(spec["path"])
    path = root / "tables/rna_effect_profiles.npz"
    with np.load(path) as z:
        arrays = {k:z[k].copy() for k in z.files}
    arrays["log2fc"][0, 0] = np.nan
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match="declared direct-target"):
        comparison.read_bundle(renew(root, "chip01"))


def test_missing_chip_is_not_a_complete_overall_result(tmp_path):
    path, cfg = config(tmp_path)
    cfg["expected_chips"].append("chip03")
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="chip coverage"):
        comparison.run(path, smoke=True)


def test_smoke_report_hashes_source_preservation_and_no_overwrite(tmp_path):
    path, cfg = config(tmp_path)
    before = {str(p):comparison.sha256(p) for root in [Path(s["path"]) for s in cfg["inputs"]]
              for p in root.rglob("*") if p.is_file()}
    result = comparison.run(path, smoke=True)
    out = Path(cfg["output"])
    assert result["scientific_result"] is False and result["counts_read"] is False
    assert comparison.verify(out)["complete"]
    assert all(comparison.sha256(p) == digest for p, digest in before.items())
    assert len(list((out/"figures").glob("*.png"))) == 3
    with pytest.raises(FileExistsError):
        comparison.run(path, smoke=True)
    (out/"unexpected.txt").write_text("fixture")
    with pytest.raises(ValueError, match="coverage"):
        comparison.verify(out)


def test_formal_comparison_requires_input_bound_review(tmp_path, monkeypatch):
    path, cfg = config(tmp_path)
    monkeypatch.setattr(comparison, "revision", lambda *a:dict(commit="a"*40, working_tree_dirty=False))
    review = tmp_path/"review.json"
    review.write_text(json.dumps(dict(compatible=True, basis="fixture", input_manifests={})))
    cfg["comparability_review"] = dict(path=str(review), sha256=comparison.sha256(review))
    path.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="missing or stale"):
        comparison.run(path)
