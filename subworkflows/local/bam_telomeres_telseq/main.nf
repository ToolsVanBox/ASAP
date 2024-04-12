//
// BAM TELOMERES TELSEQ

// Include local modules
include { TELSEQ } from '../../../modules/local/telseq/main.nf'


workflow BAM_TELOMERES_TELSEQ {
  take:
    ch_bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
  
  main:
    ch_versions = Channel.empty()
    
    TELSEQ( ch_bam_bai )
    ch_versions = ch_versions.mix( TELSEQ.out.versions )  

  emit:
    versions = ch_versions                       // channel: [ versions.yml ]

}


