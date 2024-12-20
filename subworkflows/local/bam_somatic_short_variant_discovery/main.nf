//
// BAM GERMLINE SHORT VARIANT DISCOVERY

// Include local subworkflows
include { BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2 } from '../../../subworkflows/local/bam_somatic_short_variant_discovery_mutect2/main'       

// Include nf-core modules
include { GATK4_SPLITINTERVALS } from '../../../modules/nf-core/gatk4/splitintervals/main'                                                                          

workflow BAM_SOMATIC_SHORT_VARIANT_DISCOVERY {
  take:
    ch_bam_bai_tumor  // channel: [ meta, path(bam), path(bai) ]
    ch_bam_bai_normal // channel: [ meta, path(bam), path(bai) ]
    ch_split_interval // channel: [ intervals ]
    ch_fasta // channel: [ meta, path(fasta) ]
    ch_fai // channel: [ meta, path(fai) ]
    ch_dict // channel: [ meta, path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_somatic_vcfs = Channel.empty()
    ch_somatic_tbi = Channel.empty()
    ch_somatic_f1r2 = Channel.empty()
    ch_somatic_stats = Channel.empty()
      
    for ( tool in params.bam_somatic_short_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "mutect2" ) {
        BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2( ch_bam_bai_tumor, ch_bam_bai_normal, ch_split_interval, ch_fasta, ch_fai, ch_dict )
        
        ch_versions = ch_versions.mix( BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.versions )
        ch_somatic_vcfs = ch_somatic_vcfs.mix(BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.vcf)
        ch_somatic_tbi = ch_somatic_tbi.mix(BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.tbi)
        ch_somatic_f1r2 = ch_somatic_f1r2.mix(BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.f1r2)
        ch_somatic_stats = ch_somatic_stats.mix(BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.stats)

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_somatic_vcfs       // channel: [ meta, vcf ]
    tbi = ch_somatic_tbi        // channel: [ meta, tbi ]
    f1r2 = ch_somatic_f1r2      // channel: [ meta, path(f1r2) ]
    stats = ch_somatic_stats    // channel: [ val(meta), path(stats) ]
    versions = ch_versions      // channel: [ versions.yml ]
}


