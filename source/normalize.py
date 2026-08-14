from pathlib import Path
import argparse
import logging
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def normalize_expression(adata):
    """Normalize total counts per cell and apply log1p."""
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    logger.info("Normalized total counts and applied log1p")
    return adata


def select_highly_variable_genes(adata, n_top_genes):
    """Flag highly variable genes for downstream analysis."""
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"Highly variable genes: {n_hvg}/{adata.n_vars}")
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n_top_genes", type=int, required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    adata = normalize_expression(adata)
    adata = select_highly_variable_genes(adata, args.n_top_genes)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    logger.info(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()
