//
// BAM FINGERPRINT

// Include local subworkflows
include { BAM_FINGERPRINT_PICARD } from '../../../subworkflows/local/bam_fingerprint_picard/main'                             

workflow BAM_FINGERPRINT {
  take:
    ch_bams // channel: [ meta, path(bam), path(bai) ] 
  main:
    ch_versions = Channel.empty() 
    
    def haplotype_file = file( params.genomes[params.genome].haplotype, checkIfExists: true )
    ch_haplotype = Channel.value( haplotype_file )
    
    for ( tool in params.bam_fingerprint.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Fingerprint
      if ( tool == "picard" ) {
          ch_bams_picard = ch_bams.map{ meta, bam, bai -> [ [ id:meta.run_id ], bam ] }.groupTuple()
          
          BAM_FINGERPRINT_PICARD( ch_bams_picard, [], ch_haplotype )   
          ch_versions = ch_versions.mix( BAM_FINGERPRINT_PICARD.out.versions )  
                    
          known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam fingerprint tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}