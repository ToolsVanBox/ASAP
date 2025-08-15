//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY

// Include local subworkflows
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS } from '../../../subworkflows/local/bam_germline_structural_variant_discovery_gridss/main'
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA } from '../../../subworkflows/local/bam_germline_structural_variant_discovery_manta/main'

// Include nf-core modules

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_fasta_dict // channel: [ val(meta), path(dict) ]
    
  main:
    ch_versions = Channel.empty()
    ch_gridss_vcf = Channel.empty()
    
    for ( tool in params.bam_germline_structural_variant_discovery.tool ) {
      def tool = tool.toLowerCase()      
      def known_tool = false   
      
      if ( tool == "gridss" ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS( ch_bam_bai, ch_fasta, ch_fai, ch_fasta_dict )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS.out.versions )

        ch_gridss_vcf = ch_gridss_vcf.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS.out.vcf )
      
        known_tool = true
      }

      if ( tool == "manta" ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA( ch_bam_bai )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA.out.versions )

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }

  emit:
    gridss_vcf = ch_gridss_vcf // channel: [ meta, vcf ]
    versions = ch_versions // channel: [ versions.yml ]
}


