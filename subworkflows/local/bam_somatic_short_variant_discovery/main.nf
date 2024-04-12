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
      
    for ( tool in params.bam_somatic_short_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "mutect2" ) {
        BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2( ch_bam_bai_tumor, ch_bam_bai_normal, ch_split_interval, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2.out.versions )

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_somatic_vcfs // channel: [ meta, vcf ]
    tbi = ch_somatic_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


