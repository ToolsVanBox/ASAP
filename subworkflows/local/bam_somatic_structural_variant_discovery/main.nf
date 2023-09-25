//
// BAM SOMATIC STRUCTURAL VARIANT DISCOVERY

// Include local subworkflows
include { BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX } from '../../../subworkflows/local/bam_somatic_structural_variant_discovery_gridss_purple_linx/main'
include { BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA } from '../../../subworkflows/local/bam_somatic_structural_variant_discovery_manta/main'

// Include nf-core modules

workflow BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bams  // channel: [ meta, path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai) ]
  main:
    ch_versions = Channel.empty()

    for ( tool in params.bam_somatic_structural_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false   
      
      if ( tool == "gridss_purple_linx" ) {
        BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX( ch_bams )
        ch_versions = ch_versions.mix( BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX.out.versions )
      
        known_tool = true
      }

      if ( tool == "manta" ) {
        BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA( ch_bams )
        ch_versions = ch_versions.mix( BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA.out.versions )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam somatic short variant discovery tool, this tool is not build in (yet).")
      }
    }
    
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


