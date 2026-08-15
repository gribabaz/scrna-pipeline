process CLUSTER {
    publishDir 'results/cluster', mode: 'copy'

    input:
    path norm_output

    output:
    path 'cluster_output.h5ad'

    script:
    """
    python /pipeline/source/cluster.py \
        --input ${norm_output} \
        --output cluster_output.h5ad \
        --n_pcs ${params.n_pcs} \
        --n_neighbors ${params.n_neighbors} \
        --resolution ${params.resolution}
    """
}