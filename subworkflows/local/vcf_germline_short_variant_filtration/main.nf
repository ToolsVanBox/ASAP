//
// VCF VARIANT FILTRATION

// Include local subworkflows
include { VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4 } from '../../../subworkflows/local/vcf_germline_short_variant_filtration_gatk4/main'

// Include nf-core modules

workflow VCF_GERMLINE_SHORT_VARIANT_FILTRATION {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_dict // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_vcfs = Channel.empty()
    ch_tbi = Channel.empty()

    for ( tool in params.vcf_germline_short_variant_filtration.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "gatk4" ) {
        
        VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4( ch_vcf_tbi, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4.out.versions )

        ch_vcfs = ch_vcfs.mix( VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4.out.vcf )
        ch_tbi = ch_tbi.mix( VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4.out.tbi )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf variant filtration tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_vcfs // channel: [ meta, vcf ]
    tbi = ch_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


