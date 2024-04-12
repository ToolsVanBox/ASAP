//
// BAM TELOMERES

// Include local subworkflows
include { BAM_TELOMERES_TELOMERECAT } from '../../../subworkflows/local/bam_telomeres_telomerecat/main'                             
include { BAM_TELOMERES_TELSEQ } from '../../../subworkflows/local/bam_telomeres_telseq/main'                             

workflow BAM_TELOMERES {
  take:
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ] 
  main:
    ch_versions = Channel.empty() 
    
    for ( tool in params.bam_telomeres.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Telomerecat bam2length
      if ( tool == "telomerecat" ) {

          BAM_TELOMERES_TELOMERECAT( ch_bam_bai )   
          ch_versions = ch_versions.mix( BAM_TELOMERES_TELOMERECAT.out.versions )  
                    
          known_tool = true
      }

      if ( tool == "telseq" ) {

          BAM_TELOMERES_TELSEQ( ch_bam_bai )   
          ch_versions = ch_versions.mix( BAM_TELOMERES_TELSEQ.out.versions )  
                    
          known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam telomere tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}