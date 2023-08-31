
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main.nf'


workflow PREPARE_GENOME {
    take:
        fasta
        tool
    main:

        versions = Channel.empty()

        if ( tool.toString().toUpperCase() == "BWA_MEM" ) {
            BWA_INDEX( fasta )     // If aligner is bwa-mem
            versions = versions.mix(BWA_INDEX.out.versions)
            // BWA_INDEX.out.view()
        }
    emit:
        bwa = BWA_INDEX.out.index.map{ meta, index -> [index] }.collect()      // path: bwa/*

}