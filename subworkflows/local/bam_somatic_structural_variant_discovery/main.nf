//
// BAM SOMATIC STRUCTURAL VARIANT DISCOVERY

// Include local subworkflows
include { BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA } from '../../../subworkflows/local/bam_somatic_structural_variant_discovery_manta/main'

// Include nf-core modules

workflow BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bam_bai_normal_tumor  // channel: [ meta, path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

  main:
    ch_versions = Channel.empty()

    for ( tool in params.bam_somatic_structural_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false   

      ch_manta_somatic = ch_bam_bai_normal_tumor
        .map{ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
          [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai, [], [] ]
        }
        
      if ( tool == "manta" ) {
        BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA( ch_manta_somatic, ch_fasta, ch_fai )
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


