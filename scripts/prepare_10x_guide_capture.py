"""Build a labeled RNA-count H5AD from a 10x CRISPR Guide Capture result.

This converter joins a filtered 10x feature-barcode matrix to Cell Ranger's
per-cell protospacer calls through the authoritative feature-reference table.
It deliberately keeps only singleton-guide cells by default: guide strings are
never parsed to guess a target, and multi-guide / no-guide cells are not silently
assigned to a perturbation.  The output still contains raw integer RNA counts;
run ``preprocess_h5ad.py`` separately for QC, normalization, and HVG selection.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd


REQUIRED_CALL_COLUMNS = {"cell_barcode", "num_features", "feature_call", "num_umis"}
REQUIRED_REFERENCE_COLUMNS = {"name", "target_gene_name"}


def prepare_10x_guide_capture(
    matrix_dir: Path,
    calls_path: Path,
    feature_reference_path: Path,
    output_path: Path,
    control_target: str = "Non-Targeting",
) -> Any:
    """Write singleton-guide RNA counts and resolved perturbation labels to H5AD.

    ``feature_call`` is joined exactly to ``feature_reference.name``.  A target
    label is therefore always sourced from the library reference, never inferred
    from a guide string.  Cells with zero or more than one called guide are
    excluded because they lack an unambiguous categorical perturbation label.
    """
    import scanpy as sc

    calls = pd.read_csv(calls_path)
    missing_calls = REQUIRED_CALL_COLUMNS - set(calls.columns)
    if missing_calls:
        raise KeyError(f"Missing protospacer-call columns: {sorted(missing_calls)}")
    if calls["cell_barcode"].duplicated().any():
        raise ValueError("Protospacer calls must have at most one row per cell barcode.")

    reference = pd.read_csv(feature_reference_path)
    missing_reference = REQUIRED_REFERENCE_COLUMNS - set(reference.columns)
    if missing_reference:
        raise KeyError(f"Missing feature-reference columns: {sorted(missing_reference)}")
    if reference["name"].duplicated().any():
        raise ValueError("Feature-reference guide names must be unique.")

    adata = sc.read_10x_mtx(
        matrix_dir,
        var_names="gene_symbols",
        make_unique=True,
        gex_only=False,
    )
    if "feature_types" not in adata.var:
        raise KeyError("10x feature matrix did not provide adata.var['feature_types'].")
    rna_mask = adata.var["feature_types"].astype(str).eq("Gene Expression")
    if not rna_mask.any():
        raise ValueError("No 'Gene Expression' features were found in the 10x matrix.")
    adata = adata[:, rna_mask.to_numpy()].copy()

    calls = calls.set_index("cell_barcode")
    if not calls.index.is_unique:
        raise ValueError("Protospacer calls must have at most one row per cell barcode.")
    joined_calls = calls.reindex(adata.obs_names)
    singleton_mask = joined_calls["num_features"].eq(1)
    if not singleton_mask.any():
        raise ValueError("No singleton-guide cells overlap the expression matrix.")

    adata = adata[singleton_mask.to_numpy()].copy()
    selected_calls = joined_calls.loc[adata.obs_names]
    guide_to_target = reference.set_index("name")["target_gene_name"]
    targets = selected_calls["feature_call"].map(guide_to_target)
    if targets.isna().any():
        unknown = selected_calls.loc[targets.isna(), "feature_call"].unique().tolist()
        raise ValueError(
            "Singleton guide calls are absent from feature_reference.csv: "
            f"{unknown[:10]}"
        )

    source_targets = targets.astype(str)
    adata.obs["guide_id"] = selected_calls["feature_call"].astype(str).to_numpy()
    adata.obs["source_perturbation"] = source_targets.to_numpy()
    adata.obs["perturbed_gene"] = source_targets.where(
        source_targets.ne(control_target), "ctrl"
    ).to_numpy()
    adata.obs["guide_umis"] = selected_calls["num_umis"].astype(int).to_numpy()
    adata.obs["guide_num_features"] = selected_calls["num_features"].astype(int).to_numpy()
    adata.layers["counts"] = adata.X.copy()
    adata.var["gene_name"] = adata.var_names.astype(str)

    matrix_cell_count = int(joined_calls.shape[0])
    call_count = int(joined_calls["num_features"].notna().sum())
    singleton_count = int(singleton_mask.sum())
    adata.uns["perturbseq_input"] = {
        "source_format": "10x CRISPR Guide Capture filtered feature-barcode matrix",
        "matrix_dir": str(matrix_dir),
        "protospacer_calls": str(calls_path),
        "feature_reference": str(feature_reference_path),
        "assignment_policy": "singleton-guide cells only",
        "control_target": control_target,
        "matrix_cells": matrix_cell_count,
        "cells_with_guide_calls": call_count,
        "singleton_guide_cells": singleton_count,
        "multiguide_cells": int(joined_calls["num_features"].gt(1).sum()),
        "zero_guide_matrix_cells": int(matrix_cell_count - call_count),
        "rna_feature_count": int(adata.n_vars),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Saved: {output_path}")
    print(f"Final shape: {adata.n_obs} cells x {adata.n_vars} RNA genes")
    print(f"Perturbation labels: {adata.obs['perturbed_gene'].nunique()}")
    return adata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        required=True,
        help="10x filtered_feature_bc_matrix directory.",
    )
    parser.add_argument(
        "--protospacer-calls",
        type=Path,
        required=True,
        help="Cell Ranger protospacer_calls_per_cell.csv.",
    )
    parser.add_argument(
        "--feature-reference",
        type=Path,
        required=True,
        help="Authoritative Cell Ranger feature_reference.csv.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--control-target",
        default="Non-Targeting",
        help="Reference target label that becomes canonical 'ctrl'.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    prepare_10x_guide_capture(
        matrix_dir=args.matrix_dir,
        calls_path=args.protospacer_calls,
        feature_reference_path=args.feature_reference,
        output_path=args.output,
        control_target=args.control_target,
    )
