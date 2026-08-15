process QC {
    publishDir 'results/qc', mode: 'copy'

    input:
    path raw_input

    output:
    path 'qc_output.h5ad'

    script:
    """
    python /pipeline/source/qc.py \
        --input ${raw_input} \
        --output qc_output.h5ad \
        --min_genes ${params.min_genes} \
        --max_genes ${params.max_genes} \
        --max_pct_mt ${params.max_pct_mt} \
        --min_cells ${params.min_cells} \
        --mito_prefix ${params.mito_prefix} \
        --ribo_prefix ${params.ribo_prefix.join(' ')}
    """
}