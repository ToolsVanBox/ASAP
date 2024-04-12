//
// BAM FINGERPRINT

// Include local subworkflows
include { BAM_FINGERPRINT_PICARD } from '../../../subworkflows/local/bam_fingerprint_picard/main'                             
include { BAM_FINGERPRINT_GATK } from '../../../subworkflows/local/bam_fingerprint_gatk/main' 

workflow BAM_FINGERPRINT {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ] 
    ch_fasta    // channel: [ val(meta), path(fasta) ]
    ch_fai      // channel: [ val(meta), path(fai) ]
    ch_dict     // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty() 
    
    def haplotype_file = file( params.genomes[params.genome].haplotype, checkIfExists: true )
    ch_haplotype = Channel.value( haplotype_file )
    
    for ( tool in params.bam_fingerprint.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Fingerprint
      if ( tool == "picard" ) {
          ch_bam_bai_picard = ch_bam_bai.map{ meta, bam, bai -> [ [ id:meta.run_id ], bam ] }.groupTuple()

          BAM_FINGERPRINT_PICARD( ch_bam_bai_picard, [], ch_haplotype )   
          ch_versions = ch_versions.mix( BAM_FINGERPRINT_PICARD.out.versions )  
                    
          known_tool = true
      }

      // Run GATKunifiedgenotyper
      if (tool == "gatkunifiedgenotyper") {
          //ch_bam_bai_gatk = ch_bam_bai.map{ meta, bam, bai -> [ [ id:meta.run_id ], bam, bai ] }.groupTuple()

          BAM_FINGERPRINT_GATK( ch_bam_bai, ch_fasta, ch_fai, ch_dict )

          known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam fingerprint tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}