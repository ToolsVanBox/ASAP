//
// PRE MAPPING QC

// Include nf-core subworkflows
include { BAM_QC_PICARD } from '../../../subworkflows/nf-core/bam_qc_picard/main.nf'

workflow BAM_QC_POST_MAPPING {
  take:
    ch_bams // channel: [ meta, path(bam), path(bai) ]
    
  main:
    ch_versions = Channel.empty()
    
    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( fasta.toString()+".fai", checkIfExists: true )
    def fasta_dict = file( fasta.toString()+".dict", checkIfExists: true )
    def dict = file( fasta.toString().replace(".fasta",".dict" ), checkIfExists: true )
    
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( fai )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }
    ch_dict = Channel.value( fai )
      .map{ genome_dict -> [ [ id:'dict' ], genome_dict ] }

    for ( tool in params.bam_qc_post_mapping.tool ) {
      tool = tool.toLowerCase()      
      def known_tool = false
      
      // Run Picard 
      if ( tool  == "picard" ) {
          ch_bams_picard = ch_bams.map{ meta, bam, bai -> [ meta, bam, bai, [], [] ] }

          BAM_QC_PICARD( ch_bams_picard, ch_fasta, ch_fai, ch_dict )
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
