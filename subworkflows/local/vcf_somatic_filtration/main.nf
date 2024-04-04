//
// VCF SOMATIC FILTRATION

// Include local subworkflows
include { VCF_SOMATIC_FILTRATION_SMURF } from '../../../subworkflows/local/vcf_somatic_filtration_smurf/main'

// Include nf-core modules

workflow VCF_SOMATIC_FILTRATION {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_bam_bai // channel: [ meta, bam, bai ]
    
  main:
    ch_versions = Channel.empty()

    for ( tool in params.vcf_somatic_filtration.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false 

      if ( tool == "smurf" ) {
        
        println "RUN SMURF"
        VCF_SOMATIC_FILTRATION_SMURF( ch_vcf_tbi, ch_bam_bai )
        // ch_versions = ch_versions.mix( VCF_SOMATIC_FILTRATION_SMURF.out.versions )
        
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf somatic filtration tool, this tool is not build in (yet).")
      }
    }
  emit:
    // vcf = ch_vcfs // channel: [ meta, vcf ]
    // tbi = ch_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


