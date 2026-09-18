"""Prepare a validated Perturb-seq AnnData input for the inverse-cell pipeline.

Unlike ``preprocess_adamson.py``, this entry point does not assume an Adamson
file layout.  It accepts an existing H5AD object, preserves integer RNA counts
in ``adata.layers['counts']``, records a canonical ``perturbed_gene`` column,
and writes the normalized/HVG representation expected by contrastiveVI.

The caller must explicitly identify the source perturbation column and any
control labels.  This avoids silently treating a guide ID as a target gene.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


def _matrix_values(matrix: Any) -> np.ndarray:
    """Return stored values without densifying a sparse count matrix."""
    return np.asarray(matrix.data if hasattr(matrix, "data") else matrix)


def validate_count_matrix(matrix: Any, source: str) -> None:
    """Require finite, non-negative, count-like values before model preparation."""
    values = _matrix_values(matrix)
    if not np.isfinite(values).all():
        raise ValueError(f"{source} contains non-finite values.")
    if (values < 0).any():
        raise ValueError(f"{source} contains negative values.")
    if not np.allclose(values, np.rint(values)):
        raise ValueError(
            f"{source} is not count-like. Supply the unnormalized integer-count layer, "
            "not normalized or log-transformed expression."
        )


def prepare_perturbseq_h5ad(
    input_path: Path,
    output_path: Path,
    perturbation_key: str,
    control_values: Iterable[str],
    counts_layer: str | None = None,
    guide_key: str | None = None,
    feature_type_key: str | None = None,
    gene_expression_value: str = "Gene Expression",
    gene_embedding_path: Path | None = None,
    min_genes: int = 100,
    min_gene_cells: int = 3,
    n_hvg: int = 5000,
) -> Any:
    """Normalize one Perturb-seq H5AD into the project input contract.

    ``perturbation_key`` must contain target-gene identities or already resolved
    controls.  If the source column contains guide IDs instead, resolve them
    through the authoritative guide reference before calling this function.
    """
    import scanpy as sc

    adata = sc.read_h5ad(input_path)
    if perturbation_key not in adata.obs:
        raise KeyError(f"Missing adata.obs[{perturbation_key!r}].")
    if guide_key is not None and guide_key not in adata.obs:
        raise KeyError(f"Missing adata.obs[{guide_key!r}].")
    if feature_type_key is not None and feature_type_key not in adata.var:
        raise KeyError(f"Missing adata.var[{feature_type_key!r}].")

    if feature_type_key is not None:
        rna_mask = adata.var[feature_type_key].astype(str).eq(gene_expression_value)
        if not rna_mask.any():
            raise ValueError(
                f"No features have {feature_type_key!r} == {gene_expression_value!r}."
            )
        adata = adata[:, rna_mask.to_numpy()].copy()

    source_counts = adata.X if counts_layer is None else adata.layers.get(counts_layer)
    source_name = "adata.X" if counts_layer is None else f"adata.layers[{counts_layer!r}]"
    if source_counts is None:
        raise KeyError(f"Missing {source_name}.")
    validate_count_matrix(source_counts, source_name)
    adata.layers["counts"] = source_counts.copy()

    original_labels = adata.obs[perturbation_key].astype(str)
    controls = {str(value) for value in control_values}
    adata.obs["source_perturbation"] = original_labels
    adata.obs["perturbed_gene"] = original_labels.where(~original_labels.isin(controls), "ctrl")
    if guide_key is not None:
        adata.obs["guide_id"] = adata.obs[guide_key].astype(str)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_gene_cells)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("QC filtering removed every cell or every gene.")

    if gene_embedding_path is not None:
        embeddings = pd.read_csv(gene_embedding_path, index_col=0)
        supported_targets = set(embeddings.index.astype(str)) | {"ctrl"}
        adata = adata[adata.obs["perturbed_gene"].isin(supported_targets)].copy()
        if adata.n_obs == 0:
            raise ValueError("No cells remain after restricting to GenePT-supported targets.")

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if n_hvg > 0:
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_hvg, adata.n_vars), subset=False)
        targets = set(adata.obs["perturbed_gene"].astype(str)) - {"ctrl"}
        keep = adata.var["highly_variable"].to_numpy() | adata.var_names.isin(targets)
        adata = adata[:, keep].copy()
    adata.var["gene_name"] = adata.var_names.astype(str)

    adata.uns["perturbseq_input"] = {
        "input_path": str(input_path),
        "source_counts": source_name,
        "perturbation_key": perturbation_key,
        "control_values": sorted(controls),
        "guide_key": guide_key,
        "feature_type_key": feature_type_key,
        "gene_expression_value": gene_expression_value,
        "min_genes": min_genes,
        "min_gene_cells": min_gene_cells,
        "n_hvg": n_hvg,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Saved: {output_path}")
    print(f"Final shape: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Perturbation labels: {adata.obs['perturbed_gene'].nunique()}")
    return adata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Source AnnData (.h5ad).")
    parser.add_argument("--output", type=Path, required=True, help="Prepared AnnData output.")
    parser.add_argument(
        "--perturbation-key",
        required=True,
        help="adata.obs column containing resolved perturbation target labels.",
    )
    parser.add_argument(
        "--control-value",
        action="append",
        default=[],
        help="Source perturbation value that should become canonical 'ctrl'. Repeat as needed.",
    )
    parser.add_argument(
        "--counts-layer",
        default=None,
        help="Layer containing raw integer counts; omit only when adata.X contains raw counts.",
    )
    parser.add_argument("--guide-key", default=None, help="Optional adata.obs guide-ID column.")
    parser.add_argument(
        "--feature-type-key",
        default=None,
        help="Optional adata.var column used to retain RNA features only.",
    )
    parser.add_argument("--gene-expression-value", default="Gene Expression")
    parser.add_argument("--gene-embeddings", type=Path, default=None)
    parser.add_argument("--min-genes", type=int, default=100)
    parser.add_argument("--min-gene-cells", type=int, default=3)
    parser.add_argument("--n-hvg", type=int, default=5000, help="Use 0 to skip HVG selection.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    prepare_perturbseq_h5ad(
        input_path=args.input,
        output_path=args.output,
        perturbation_key=args.perturbation_key,
        control_values=args.control_value,
        counts_layer=args.counts_layer,
        guide_key=args.guide_key,
        feature_type_key=args.feature_type_key,
        gene_expression_value=args.gene_expression_value,
        gene_embedding_path=args.gene_embeddings,
        min_genes=args.min_genes,
        min_gene_cells=args.min_gene_cells,
        n_hvg=args.n_hvg,
    )
