import argparse
import logging
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def run_pca(adata, n_pcs):
    """Compute PCA using highly variable genes."""
    sc.pp.pca(adata, n_comps=n_pcs)
    logger.info(f"PCA computed: {n_pcs} components")
    return adata


def compute_neighbors(adata, n_neighbors):
    """Build the KNN neighbors graph on PCA components."""
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    logger.info(f"Neighbors graph computed: n_neighbors={n_neighbors}")
    return adata


def run_leiden_clustering(adata, resolution):
    """Cluster cells using the Leiden algorithm."""
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
    n_clusters = adata.obs["leiden"].nunique()
    logger.info(f"Leiden clustering done: {n_clusters} clusters")
    return adata


def run_umap(adata):
    """Compute UMAP coordinates for visualization."""
    sc.tl.umap(adata)
    logger.info("UMAP computed")
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n_pcs", type=int, required=True)
    parser.add_argument("--n_neighbors", type=int, required=True)
    parser.add_argument("--resolution", type=float, required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    adata = run_pca(adata, args.n_pcs)
    adata = compute_neighbors(adata, args.n_neighbors)
    adata = run_leiden_clustering(adata, args.resolution)
    adata = run_umap(adata)

    adata.write_h5ad(args.output)
    logger.info(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()