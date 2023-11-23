//
// BAM TELOMERES TELOMERECAT

// Include local modules
include { TELOMERECAT_BAM2LENGTH } from '../../../modules/local/telomerecat/bam2length/main'                             


workflow BAM_TELOMERES_TELOMERECAT {
  take:
    ch_bam_bai  // channel: [ val(meta), path(bam) ]
  
  main:
    ch_versions = Channel.empty()
    
    TELOMERECAT_BAM2LENGTH( ch_bam_bai )
    ch_versions = ch_versions.mix( TELOMERECAT_BAM2LENGTH.out.versions )  

  emit:
    versions = ch_versions                       // channel: [ versions.yml ]

}


