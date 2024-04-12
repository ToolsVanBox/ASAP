//
// Create genome index files
//

// Include nf-core modules
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main.nf'
include { BWAMEM2_INDEX } from '../../../modules/nf-core/bwamem2/index/main'                                   

workflow PREPARE_GENOME {
    take:
        fasta // channel: [ meta, fasta ]
        tool // val(tool)
    main:
        versions = Channel.empty()

        if ( tool.toString().toLowerCase() == "bwamem1" ) {
            BWA_INDEX( fasta )     // If aligner is bwa-mem
            versions = versions.mix(BWA_INDEX.out.versions)
            index = BWA_INDEX.out.index.map{ meta, index -> [index] }.collect()      // path: bwa/*
        } 
        else if ( tool.toString().toLowerCase() == "bwamem2" ) {
            BWAMEM2_INDEX( fasta )     // If aligner is bwa-mem2
            versions = versions.mix(BWAMEM2_INDEX.out.versions)
            index = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()      // path: bwamem2/*
        }
     emit:
        index = index // channel: path: bwa/*
        versions = versions // channel: [ versions.yml ]

}