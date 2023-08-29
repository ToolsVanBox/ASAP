
nextflow.enable.dsl=2

include { WAP } from './workflows/wap.nf'

workflow {
    WAP()
}