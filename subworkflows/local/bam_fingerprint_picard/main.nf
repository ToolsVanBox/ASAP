//
// BAM FINGERPRINT PICARD

include { PICARD_CROSSCHECKFINGERPRINTS } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'                             


workflow BAM_FINGERPRINT_PICARD {
  take:
    ch_bams   // channel: [ val(meta), path(bam) ]
    ch_input2 // channel: [ path(input2) ]
    ch_haplotype // channel: [ path(haplotype) ]    
  
  main:

    ch_versions = Channel.empty()
    println "RUN PICARD CROSS CHECK FINGERPRINTS"
    
    PICARD_CROSSCHECKFINGERPRINTS( ch_bams, ch_input2, ch_haplotype )

  emit:
    versions = ch_versions                       // channel: [ versions.yml ]

}


