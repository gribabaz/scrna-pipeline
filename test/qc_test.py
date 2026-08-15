import numpy as np
import pandas as pd
import anndata as ad
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "source"))

from qc import calculate_qc_metrics, filter_cells, filter_genes

@pytest.fixture
def small_adata():
    """Minimal synthetic AnnData for QC threshold testing"""
    gene_names = ["MT-CO1", "MT-ND1", "RPS4", "GENE_A", "GENE_B", "GENE_C"]

    X = np.array([
        [10,  5,  3, 100, 100, 100],
        [50, 40,  1,   5,   0,   0],
        [ 0,  0,  0,   1,   0,   0],
        [ 5,  5,  5,  50,  50,  50],
        [ 5,  5,  5,  50,  50,  50],
    ])

    adata = ad.AnnData(
        X=X.astype(float),
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(5)]),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


def test_calculate_qc_metrics(small_adata):
    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix=["RPS"])
    assert adata.obs["pct_counts_mt"].iloc[1] > 0
    assert adata.obs["pct_counts_mt"].iloc[2] == 0


def test_filter_cells(small_adata):
    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix=["RPS"])
    adata = filter_cells(adata, min_genes=2, max_genes=10, max_pct_mt=20.0)

    assert adata.n_obs == 3
    assert adata.obs_names.tolist() == ["cell_0", "cell_3", "cell_4"]


def test_filter_genes(small_adata, min_cells=3):
    adata = filter_genes(small_adata)
    assert "GENE_B" in adata.var_names
    assert adata.n_vars <= small_adata.n_vars
