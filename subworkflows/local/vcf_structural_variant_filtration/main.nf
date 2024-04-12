//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY WITH GRIDSS PURPLE LINX

// Include local subworkflows
include { VCF_STRUCTURAL_VARIANT_FILTRATION_GRIDSS_PURPLE_LINX } from '../../../subworkflows/local/vcf_structural_variant_filtration_gridss_purple_linx/main'

// Include nf-core modules

workflow VCF_STRUCTURAL_VARIANT_FILTRATION {
  take:
    ch_gridss_vcf  // channel: [ meta, path(vcf) ]
    ch_bam_bai_normal // channel: [ meta, path(bam), path(bai) ]
    ch_bam_bai_tumor // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_fasta_dict // channel: [ val(meta), path(dict) ]
    
  main:
    ch_versions = Channel.empty()
    
    for ( tool in params.vcf_structural_variant_filtration.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "gridss_purple_linx" ) {
        
        VCF_STRUCTURAL_VARIANT_FILTRATION_GRIDSS_PURPLE_LINX( ch_gridss_vcf, ch_bam_bai_normal, ch_bam_bai_tumor, ch_fasta, ch_fai, ch_fasta_dict )
        ch_versions = ch_versions.mix( VCF_STRUCTURAL_VARIANT_FILTRATION_GRIDSS_PURPLE_LINX.out.versions )
        
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf structural variant filtration tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


