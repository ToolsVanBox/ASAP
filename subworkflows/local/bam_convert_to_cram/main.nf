//
// BAM CONVERT TO CRAM

// Include nf-core modules
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'                             

workflow BAM_CONVERT_TO_CRAM {
  take:
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ] 
    ch_fai // channel: [ val(meta), path(fai) ] 

  main:
    ch_versions = Channel.empty() 
    ch_crams = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    fai = ch_fai.map{ meta, fai -> [ fai ] }

    for ( tool in params.bam_convert_to_cram.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Fingerprint
      if ( tool == "samtools" ) {
        SAMTOOLS_CONVERT( ch_bam_bai, fasta, fai )
        ch_versions = ch_versions.mix( SAMTOOLS_CONVERT.out.versions )  

        ch_crams = SAMTOOLS_CONVERT.out.alignment_index
                    
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam convert to cram tool, this tool is not build in (yet).")
      }
    }
  emit:
    cram_crai = ch_crams // channel: [ meta, cram, crai ]
    versions = ch_versions // channel: [ versions.yml ]
}