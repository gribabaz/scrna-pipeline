import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "source"))

from normalize import normalize_expression, select_highly_variable_genes


@pytest.fixture
def small_adata():
    """AnnData for normalization unit tests"""
    rng = np.random.default_rng(0)
    X_raw = rng.poisson(3, size=(50, 200)).astype(float)
    return ad.AnnData(
        X=X_raw,
        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(50)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(200)]),
    )


def test_normalize_expression(small_adata):
    X_raw = small_adata.X.copy()

    adata = normalize_expression(small_adata)

    assert not np.allclose(adata.X, X_raw)
    assert adata.X.min() >= 0


def test_select_highly_variable_genes(small_adata):
    adata = normalize_expression(small_adata)
    adata = select_highly_variable_genes(adata, n_top_genes=20)

    assert adata.var["highly_variable"].sum() == 20