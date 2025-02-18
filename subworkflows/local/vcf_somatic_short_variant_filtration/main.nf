//
// VCF SOMATIC SHORT VARIANT FILTRATION
//

// Include local subworkflows
include { VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2 } from '../../../subworkflows/local/vcf_somatic_short_variant_filtration_mutect2/main'

// Include nf-core modules

workflow VCF_SOMATIC_SHORT_VARIANT_FILTRATION {
  take:
    ch_somatic_vcf          // channel: [ val(meta), path(vcf) ]
    ch_somatic_tbi          // channel: [ val(meta), path(tbi) ]
    ch_f1r2                 // channel: [ val(meta), path(f1r2) ]
    ch_somatic_stats        // channel: [ val(meta), path(stats) ]
    ch_bam_bai_tumor        // channel: [ val(meta), path(bam), path(bai) ]
    ch_bam_bai_normal       // channel: [ val(meta), path(bam), path(bai) ]
    ch_split_intervals      // channel: [ intervals ]
    ch_fasta                // channel: [ meta, path(fasta) ]
    ch_fai                  // channel: [ meta, path(fai) ]
    ch_dict                 // channel: [ meta, path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_vcfs = Channel.empty()
    ch_tbi = Channel.empty()
    ch_stats = Channel.empty()

    for ( tool in params.vcf_somatic_short_variant_filtration.tool ) {
      def tool = tool.toLowerCase()      
      def known_tool = false    

      if ( tool == "mutect2" ) {
        
        VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2( ch_somatic_vcf, ch_somatic_tbi, ch_f1r2, ch_somatic_stats, ch_bam_bai_tumor, ch_bam_bai_normal, ch_split_intervals, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2.out.versions )

        ch_vcfs = ch_vcfs.mix( VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2.out.filtered_vcf )
        ch_tbi = ch_tbi.mix( VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2.out.filtered_tbi )
        ch_stats = ch_stats.mix( VCF_SOMATIC_SHORT_VARIANT_FILTRATION_MUTECT2.out.filtered_stats )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf variant filtration tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_vcfs                   // channel: [ meta, vcf ]
    tbi = ch_tbi                    // channel: [ meta, tbi ]
    stats = ch_stats                // channel: [ meta, stats ]
    versions = ch_versions          // channel: [ versions.yml ]
}