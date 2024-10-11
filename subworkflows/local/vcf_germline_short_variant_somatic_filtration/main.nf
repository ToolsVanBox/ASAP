//
// VCF SOMATIC FILTRATION

// Include local subworkflows
include { VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF } from '../../../subworkflows/local/vcf_germline_short_variant_somatic_filtration_smurf/main'

// Include nf-core modules

workflow VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_bam_bai // channel: [ meta, bam, bai ]
    
  main:
    ch_versions = Channel.empty()

    for ( tool in params.vcf_germline_short_variant_somatic_filtration.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false 

      if ( tool == "smurf" ) {
        
        VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF( ch_vcf_tbi, ch_bam_bai )
        ch_filtered_vcfs = VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF.out.filtered_vcf
        ch_filtered_tbi = VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF.out.filtered_tbi
        ch_vcfs = VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF.out.vcf
        ch_tbi = VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF.out.tbi
        ch_versions = ch_versions.mix( VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF.out.versions )
        
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf somatic filtration tool, this tool is not build in (yet).")
      }
    }
  emit:
    filtered_vcf = ch_filtered_vcfs   // channel: [ meta, vcf ]
    filtered_tbi = ch_filtered_tbi    // channel: [ meta, tbi ]
    vcf = ch_vcfs                     // channel: [ meta, vcf ]
    tbi = ch_tbi                      // channel: [ meta, tbi ] 
    versions = ch_versions            // channel: [ versions.yml ]
}


