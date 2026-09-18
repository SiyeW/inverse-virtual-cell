"""Integration checks for 10x CRISPR Guide Capture conversion."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import read_h5ad
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_10x_guide_capture import prepare_10x_guide_capture


def _write_gzip_text(path: Path, content: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(content)


def test_prepare_10x_guide_capture_uses_reference_and_singletons(tmp_path: Path) -> None:
    matrix_dir = tmp_path / "filtered_feature_bc_matrix"
    matrix_dir.mkdir()
    matrix = csr_matrix(
        np.array(
            [
                [3, 0, 5, 2],
                [0, 4, 1, 3],
                [7, 8, 9, 0],
            ],
            dtype=np.int64,
        )
    )
    temporary_matrix = tmp_path / "matrix.mtx"
    mmwrite(temporary_matrix, matrix)
    with temporary_matrix.open("rb") as source, gzip.open(matrix_dir / "matrix.mtx.gz", "wb") as sink:
        sink.write(source.read())
    _write_gzip_text(matrix_dir / "barcodes.tsv.gz", "cell_a\ncell_b\ncell_c\ncell_d\n")
    _write_gzip_text(
        matrix_dir / "features.tsv.gz",
        "gene_a\tGENE_A\tGene Expression\n"
        "gene_b\tGENE_B\tGene Expression\n"
        "guide_1\tGUIDE_1\tCRISPR Guide Capture\n",
    )

    calls_path = tmp_path / "protospacer_calls_per_cell.csv"
    pd.DataFrame(
        {
            "cell_barcode": ["cell_a", "cell_b", "cell_c"],
            "num_features": [1, 1, 2],
            "feature_call": ["GENE_A-1", "Non-Targeting-1", "GENE_A-1|Non-Targeting-1"],
            "num_umis": [5, 6, 7],
        }
    ).to_csv(calls_path, index=False)
    reference_path = tmp_path / "feature_reference.csv"
    pd.DataFrame(
        {
            "name": ["GENE_A-1", "Non-Targeting-1"],
            "target_gene_name": ["GENE_A", "Non-Targeting"],
        }
    ).to_csv(reference_path, index=False)

    output_path = tmp_path / "labeled.h5ad"
    prepared = prepare_10x_guide_capture(
        matrix_dir=matrix_dir,
        calls_path=calls_path,
        feature_reference_path=reference_path,
        output_path=output_path,
    )
    reread = read_h5ad(output_path)

    assert prepared.obs_names.tolist() == ["cell_a", "cell_b"]
    assert prepared.var_names.tolist() == ["GENE_A", "GENE_B"]
    assert prepared.obs["guide_id"].tolist() == ["GENE_A-1", "Non-Targeting-1"]
    assert prepared.obs["perturbed_gene"].tolist() == ["GENE_A", "ctrl"]
    assert prepared.layers["counts"].astype(int).toarray().tolist() == [[3, 0], [0, 4]]
    assert reread.uns["perturbseq_input"]["multiguide_cells"] == 1
    assert reread.uns["perturbseq_input"]["zero_guide_matrix_cells"] == 1
