
nextflow.enable.dsl=2

include { ASAP } from './workflows/asap.nf'

workflow {
    ASAP()
}