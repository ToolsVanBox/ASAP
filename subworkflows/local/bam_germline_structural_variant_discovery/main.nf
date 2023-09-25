//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY

// Include local subworkflows
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX } from '../../../subworkflows/local/bam_germline_structural_variant_discovery_gridss_purple_linx/main'
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA } from '../../../subworkflows/local/bam_germline_structural_variant_discovery_manta/main'

// Include nf-core modules

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bams  // channel: [ meta, path(bam), path(bai) ]
  main:
    ch_versions = Channel.empty()

    for ( tool in params.bam_germline_structural_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false   
      
      if ( tool == "gridss_purple_linx" ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX( ch_bams )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX.out.versions )
      
        known_tool = true
      }

      if ( tool == "manta" ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA( ch_bams )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA.out.versions )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }

  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


