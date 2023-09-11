//
// BAM FINGERPRINT PICARD

// Include nf-core modules
include { PICARD_CROSSCHECKFINGERPRINTS } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'                             


workflow BAM_FINGERPRINT_PICARD {
  take:
    ch_bams   // channel: [ val(meta), path(bam) ]
    ch_input2 // channel: [ path(input2) ]
    ch_haplotype // channel: [ path(haplotype) ]    
  
  main:
    ch_versions = Channel.empty()
    
    PICARD_CROSSCHECKFINGERPRINTS( ch_bams, ch_input2, ch_haplotype )
    ch_versions = ch_versions.mix( PICARD_CROSSCHECKFINGERPRINTS.out.versions )  

  emit:
    versions = ch_versions                       // channel: [ versions.yml ]

}


