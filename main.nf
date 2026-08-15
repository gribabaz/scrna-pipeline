nextflow.enable.dsl = 2

include { QC } from './modules/qc.nf'
include { NORMALIZE } from './modules/normalize.nf'
include { CLUSTER } from './modules/cluster.nf'

workflow {
    input_ch = Channel.fromPath(params.input)

    qc_ch   = QC(input_ch)
    norm_ch = NORMALIZE(qc_ch)
    CLUSTER(norm_ch)
}