process NORMALIZE {
    publishDir 'results/normalize', mode: 'copy'

    input:
    path qc_output

    output:
    path 'norm_output.h5ad'

    script:
    """
    python /pipeline/source/normalize.py \
        --input ${qc_output} \
        --output norm_output.h5ad \
        --n_top_genes ${params.n_top_genes}
    """
}