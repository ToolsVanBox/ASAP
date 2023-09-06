//
// BAM FINGERPRINT

include { BAM_FINGERPRINT_PICARD } from '../../local/bam_fingerprint_picard/main'                             


workflow BAM_FINGERPRINT {
  take:
    input_bams 
  main:
    
    // Convert input tools to UpperCase
    bam_fingerprint_tools = params.bam_fingerprint.tool.toString().toUpperCase()

    ch_bams = input_bams.map{ meta, bam, bai -> [ meta + [ step:"bam_fingerprint" ], bam ] } 
    
    def haplotype_file = file( params.genomes[params.genome].haplotype, checkIfExists: true )
    
    ch_haplotype = Channel.value( haplotype_file )
    ch_versions = Channel.empty()

    if ( bam_fingerprint_tools == "PICARD" | bam_fingerprint_tools.contains("PICARD") ) {
      tool = "picard"
      println "${tool}"
      ch_bams = ch_bams
        .map{ meta, bam ->
          [ [], bam ]
        }
        .groupTuple()
        
      // ch_input2 = Channel.value( null )
      
      BAM_FINGERPRINT_PICARD( ch_bams, [], ch_haplotype )   
      ch_versions = ch_versions.mix(BAM_FINGERPRINT_PICARD.out.versions)

    }

}


