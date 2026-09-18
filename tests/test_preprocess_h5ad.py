"""Integration checks for the generic H5AD preprocessing entry point."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preprocess_h5ad import prepare_perturbseq_h5ad


def test_prepare_h5ad_preserves_counts_and_labels(tmp_path: Path) -> None:
    source = AnnData(
        X=np.array(
            [
                [2, 0, 9],
                [0, 3, 8],
                [1, 1, 7],
                [4, 0, 6],
            ],
            dtype=np.int64,
        ),
        obs=pd.DataFrame(
            {
                "target": ["GENE1", "NTC", "GENE1", "GENE2"],
                "guide": ["GENE1-1", "NTC-1", "GENE1-2", "GENE2-1"],
            },
            index=["cell_a", "cell_b", "cell_c", "cell_d"],
        ),
        var=pd.DataFrame(
            {"feature_type": ["Gene Expression", "Gene Expression", "CRISPR Guide Capture"]},
            index=["GENE1", "GENE2", "GUIDE_FEATURE"],
        ),
    )
    input_path = tmp_path / "source.h5ad"
    output_path = tmp_path / "prepared.h5ad"
    source.write_h5ad(input_path)

    prepared = prepare_perturbseq_h5ad(
        input_path=input_path,
        output_path=output_path,
        perturbation_key="target",
        control_values=["NTC"],
        guide_key="guide",
        feature_type_key="feature_type",
        min_genes=1,
        min_gene_cells=1,
        n_hvg=0,
    )

    assert output_path.exists()
    assert prepared.var_names.tolist() == ["GENE1", "GENE2"]
    assert prepared.layers["counts"].astype(int).tolist() == [[2, 0], [0, 3], [1, 1], [4, 0]]
    assert prepared.obs["perturbed_gene"].tolist() == ["GENE1", "ctrl", "GENE1", "GENE2"]
    assert prepared.obs["guide_id"].tolist() == ["GENE1-1", "NTC-1", "GENE1-2", "GENE2-1"]
    assert prepared.uns["perturbseq_input"]["source_counts"] == "adata.X"
