//
// PRE MAPPING QC

// Include nf-core subworkflows
include { BAM_QC_PICARD } from '../../../subworkflows/nf-core/bam_qc_picard/main.nf'

workflow BAM_QC_POST_MAPPING {
  take:
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    
  main:
    ch_versions = Channel.empty()

    for ( tool in params.bam_qc_post_mapping.tool ) {
      tool = tool.toLowerCase()      
      def known_tool = false
      
      // Run Picard 
      if ( tool  == "picard" ) {
          ch_bam_bai_picard = ch_bam_bai.map{ meta, bam, bai -> [ meta, bam, bai, [], [] ] }

          BAM_QC_PICARD( ch_bam_bai_picard, ch_fasta, ch_fai, ch_dict )
          ch_versions = ch_versions.mix( BAM_QC_PICARD.out.versions.first() )

          known_tool = true
      }      
      
      if ( ! known_tool ) {
        println ("WARNING: Skip ${tool}, because it's not known as a QC post mapping tool, this tool is not build in (yet).")
      }
  }

  emit:
    versions = ch_versions // channel: [ versions.yml ]
        
}
