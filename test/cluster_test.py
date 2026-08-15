import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "source"))

from cluster import run_pca, compute_neighbors, run_leiden_clustering, run_umap


@pytest.fixture
def normalized_adata():
    rng = np.random.default_rng(0)
    X = rng.poisson(2, size=(200, 100)).astype(float)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(200)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(100)]),
    )
    adata.var["highly_variable"] = True
    return adata


def test_run_pca(normalized_adata):
    adata = run_pca(normalized_adata, n_pcs=10)
    assert adata.obsm["X_pca"].shape[1] == 10


def test_compute_neighbors(normalized_adata):
    adata = run_pca(normalized_adata, n_pcs=10)
    adata = compute_neighbors(adata, n_neighbors=10)
    assert "neighbors" in adata.uns


def test_run_leiden_clustering(normalized_adata):
    adata = run_pca(normalized_adata, n_pcs=10)
    adata = compute_neighbors(adata, n_neighbors=10)
    adata = run_leiden_clustering(adata, resolution=1.0)
    assert "leiden" in adata.obs.columns
    assert adata.obs["leiden"].nunique() > 0


def test_run_umap(normalized_adata):
    adata = run_pca(normalized_adata, n_pcs=10)
    adata = compute_neighbors(adata, n_neighbors=10)
    adata = run_umap(adata)
    assert adata.obsm["X_umap"].shape[1] == 2