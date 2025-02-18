//
// BAM TUMOR ONLY STRUCTURAL VARIANT DISCOVERY

// Include local subworkflows
include { BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_MANTA } from '../../../subworkflows/local/bam_tumoronly_structural_variant_discovery_manta/main'

// Include nf-core modules

workflow BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bam_bai_tumor  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

  main:
    ch_versions = Channel.empty()

    for ( tool in params.bam_tumoronly_structural_variant_discovery.tool ) {
      def tool = tool.toLowerCase()      
      def known_tool = false   
      
      if ( tool == "manta" ) {
        BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_MANTA( ch_bam_bai_tumor, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_MANTA.out.versions )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam tumor only short variant discovery tool, this tool is not build in (yet).")
      }
    }
    
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


