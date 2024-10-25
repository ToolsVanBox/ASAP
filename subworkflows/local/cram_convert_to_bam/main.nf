//
// CRAM CONVERT TO BAM

// Include nf-core modules
include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'                             

workflow CRAM_CONVERT_TO_BAM {
  take:
    ch_cram_crai // channel: [ meta, path(cram), path(crai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ] 
    ch_fai // channel: [ val(meta), path(fai) ] 

  main:
    ch_versions = Channel.empty() 
    ch_crams = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    fai = ch_fai.map{ meta, fai -> [ fai ] }

    for ( tool in params.cram_convert_to_bam.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Fingerprint
      if ( tool == "samtools" ) {
        SAMTOOLS_CONVERT( ch_cram_crai, fasta, fai )
        ch_versions = ch_versions.mix( SAMTOOLS_CONVERT.out.versions )  

        ch_bams = SAMTOOLS_CONVERT.out.alignment_index
                    
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a cram convert to bam tool, this tool is not build in (yet).")
      }
    }
  emit:
    bam_bai = ch_bams // channel: [ meta, cram, crai ]
    versions = ch_versions // channel: [ versions.yml ]
}