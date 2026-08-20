from pathlib import Path
import argparse
import logging
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def load_data(path):
    """Load 10x data from a directory (mtx) or a single .h5 file."""
    if Path(path).is_dir():
        adata = sc.read_10x_mtx(path, var_names="gene_symbols")
    else:
        adata = sc.read_10x_h5(path)
        
    adata.var_names_make_unique()
    logger.info(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes") 
    return adata


def calculate_qc_metrics(adata, mito_prefix, ribo_prefix):
    """Compute  mitochondrial and ribosomal per-cell QC metrics"""
    var_names_upper = adata.var_names.str.upper()
    adata.var["mt"] = var_names_upper.str.startswith(mito_prefix.upper())
    adata.var["ribo"] = var_names_upper.str.startswith(tuple(ribo_prefix))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True,
    )

    logger.info(
        f"median n_genes={adata.obs['n_genes_by_counts'].median():.0f}, "
        f"median pct_mito={adata.obs['pct_counts_mt'].median():.2f}%"
    )
    return adata


def filter_cells(adata, min_genes, max_genes, max_pct_mt):
    """Filter cells by gene count and mito percentage thresholds."""
    logger.info(f"Filtering {adata.n_obs} cells")
    mask = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    )
    adata = adata[mask].copy()
    logger.info(f"Cells remaining: {adata.n_obs}")
    return adata


def filter_genes(adata, min_cells):
    """Filter out unexpressed genes."""
    logger.info(f"Filtering {adata.n_vars} genes")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    logger.info(f"Genes remaining: {adata.n_vars}")
    return adata


def detect_doublets(adata):
    """Filter cell doublets using Scrublet."""
    sc.pp.scrublet(adata)
    n_doublets = adata.obs["predicted_doublet"].sum()
    logger.info(f"Doublets filtered: {n_doublets}/{adata.n_obs}")
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    return adata


#==========main=========
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min_genes", type=int, required=True)
    parser.add_argument("--max_genes", type=int, required=True)
    parser.add_argument("--max_pct_mt", type=float, required=True)
    parser.add_argument("--min_cells", type=int, required=True)
    parser.add_argument("--mito_prefix", required=True)
    parser.add_argument("--ribo_prefix", nargs="+", required=True)
    args = parser.parse_args()

    # load and annotate raw data
    adata = load_data(args.input)
    adata = calculate_qc_metrics(adata, args.mito_prefix, args.ribo_prefix)

    # filter cells, doublets, genes
    adata = filter_cells(adata, args.min_genes, args.max_genes, args.max_pct_mt)
    adata = detect_doublets(adata)
    adata = filter_genes(adata, args.min_cells)

    # save output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    logger.info(f"Finished. Saved to: {args.output}")


if __name__ == "__main__":
    main()