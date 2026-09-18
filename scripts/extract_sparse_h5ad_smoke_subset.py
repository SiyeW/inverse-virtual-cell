"""Extract a small raw-count H5AD subset without loading a dense source ``X``.

KOLF Strong Perturbations stores its processed expression matrix ``X`` densely.
This utility instead reads only selected rows and columns from an AnnData CSC
``layers['counts']`` matrix through HDF5, so a local smoke test can remain small
even when the complete source H5AD is tens of gigabytes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csc_matrix


def _decode(values: Iterable[Any]) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _categorical_values(group: Any, indices: np.ndarray | None = None) -> np.ndarray:
    categories = np.asarray(_decode(group["categories"][:]), dtype=object)
    codes = group["codes"][:] if indices is None else group["codes"][indices]
    if (codes < 0).any():
        raise ValueError("Required categorical metadata contains missing values.")
    return categories[codes]


def _index_values(parent: Any, indices: np.ndarray | None = None) -> list[str]:
    """Read an AnnData dataframe index stored as strings or as a categorical."""
    index = parent["_index"]
    if hasattr(index, "shape"):
        values = index[:] if indices is None else index[indices]
        return _decode(values)
    if "values" in index:
        values = index["values"][:] if indices is None else index["values"][indices]
        return _decode(values)
    values = _categorical_values(index, indices)
    return values.astype(str).tolist()


def _read_csc_subset(counts: Any, row_indices: np.ndarray, column_indices: np.ndarray) -> csc_matrix:
    """Read selected rows/columns from an H5AD CSC matrix without densifying it."""
    shape = tuple(int(value) for value in counts.attrs["shape"])
    encoding_type = counts.attrs.get("encoding-type", "")
    if isinstance(encoding_type, bytes):
        encoding_type = encoding_type.decode()
    if encoding_type != "csc_matrix":
        raise ValueError("Expected layers['counts'] to be encoded as a CSC matrix.")
    row_map = np.full(shape[0], -1, dtype=np.int32)
    row_map[row_indices] = np.arange(len(row_indices), dtype=np.int32)
    output_data: list[np.ndarray] = []
    output_rows: list[np.ndarray] = []
    output_indptr = [0]
    source_indptr = counts["indptr"]
    source_rows = counts["indices"]
    source_data = counts["data"]
    for column in column_indices:
        start, end = int(source_indptr[column]), int(source_indptr[column + 1])
        mapped_rows = row_map[source_rows[start:end]]
        keep = mapped_rows >= 0
        output_rows.append(mapped_rows[keep])
        output_data.append(source_data[start:end][keep])
        output_indptr.append(output_indptr[-1] + int(keep.sum()))
    data = np.concatenate(output_data) if output_data else np.array([], dtype=np.float32)
    rows = np.concatenate(output_rows) if output_rows else np.array([], dtype=np.int32)
    return csc_matrix((data, rows, np.asarray(output_indptr)), shape=(len(row_indices), len(column_indices))).tocsr()


def extract_smoke_subset(
    input_path: Path,
    output_path: Path,
    targets: list[str],
    cells_per_label: int = 64,
    n_hvg: int = 500,
    control_label: str = "NTC",
    seed: int = 42,
) -> AnnData:
    """Create a balanced raw-count KOLF-like smoke subset from a backed H5AD."""
    import h5py

    labels_to_keep = [control_label, *targets]
    if len(set(labels_to_keep)) != len(labels_to_keep):
        raise ValueError("Control and target labels must be distinct.")
    with h5py.File(input_path, "r") as source:
        required_obs = {"gene_target", "perturbation", "batch"}
        missing_obs = required_obs - set(source["obs"].keys())
        if missing_obs:
            raise KeyError(f"Missing required obs fields: {sorted(missing_obs)}")
        labels = _categorical_values(source["obs"]["gene_target"])
        rng = np.random.default_rng(seed)
        selected_rows: list[int] = []
        for label in labels_to_keep:
            candidates = np.flatnonzero(labels == label)
            if len(candidates) < cells_per_label:
                raise ValueError(f"{label!r} has {len(candidates)} cells, needs {cells_per_label}.")
            selected_rows.extend(rng.choice(candidates, size=cells_per_label, replace=False).tolist())
        selected_rows_array = np.sort(np.asarray(selected_rows, dtype=np.int64))

        var_names = np.asarray(_index_values(source["var"]), dtype=object)
        hvg_mask = source["var"]["highly_variable"][:].astype(bool)
        selected_columns = np.flatnonzero(hvg_mask)[:n_hvg].tolist()
        for target in targets:
            matches = np.flatnonzero(var_names == target)
            if len(matches):
                selected_columns.append(int(matches[0]))
        selected_columns_array = np.asarray(sorted(set(selected_columns)), dtype=np.int64)
        if not len(selected_columns_array):
            raise ValueError("No genes selected for smoke subset.")

        counts = _read_csc_subset(source["layers"]["counts"], selected_rows_array, selected_columns_array)
        obs = pd.DataFrame(
            {
                "source_perturbation": _categorical_values(source["obs"]["gene_target"], selected_rows_array),
                "guide_id": _categorical_values(source["obs"]["perturbation"], selected_rows_array),
                "batch": _categorical_values(source["obs"]["batch"], selected_rows_array),
            },
            index=_index_values(source["obs"], selected_rows_array),
        )
        obs["perturbed_gene"] = obs["source_perturbation"].replace({control_label: "ctrl"})
        var = pd.DataFrame(
            {"highly_variable_source": hvg_mask[selected_columns_array]},
            index=var_names[selected_columns_array],
        )
        var["gene_name"] = var.index.astype(str)

    subset = AnnData(X=counts.copy(), obs=obs, var=var)
    subset.layers["counts"] = counts.copy()
    subset.uns["smoke_subset"] = {
        "source_path": str(input_path),
        "source_counts_layer": "counts",
        "selection": "uniform random within requested gene_target labels",
        "control_label": control_label,
        "targets": targets,
        "cells_per_label": cells_per_label,
        "n_hvg_requested": n_hvg,
        "n_genes_selected": int(subset.n_vars),
        "seed": seed,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(output_path)
    print(f"Saved: {output_path}")
    print(f"Final shape: {subset.n_obs} cells x {subset.n_vars} genes")
    return subset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", action="append", required=True)
    parser.add_argument("--cells-per-label", type=int, default=64)
    parser.add_argument("--n-hvg", type=int, default=500)
    parser.add_argument("--control-label", default="NTC")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    extract_smoke_subset(
        input_path=args.input,
        output_path=args.output,
        targets=args.target,
        cells_per_label=args.cells_per_label,
        n_hvg=args.n_hvg,
        control_label=args.control_label,
        seed=args.seed,
    )
