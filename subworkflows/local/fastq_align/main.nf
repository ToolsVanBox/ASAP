//
// FASTQC

include { FASTQ_ALIGN_BWA } from '../../nf-core/fastq_align_bwa/main.nf'
include { PREPARE_GENOME } from '../prepare_genome/main.nf'

workflow FASTQ_ALIGN {
  take:
    input_fastqs
  main:

      // Convert input tools to UpperCase
    fastqc_align_tools = params.fastq_align.tool.toString().toUpperCase()
    
    if ( fastqc_align_tools == "BWA_MEM" | fastqc_align_tools.contains("BWA_MEM") ) {
        // tool = "BWA_MEM"
        genome_fasta = params.genomes[params.genome].fasta
        
        // PREPARE_GENOME( ch_fasta, tool )
        
        // PREPARE_GENOME.out.bwa.view() 

        bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
          .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
          .view()

        bam_sort = true
        ch_fasta = Channel.fromPath( genome_fasta )
          .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }

        FASTQ_ALIGN_BWA( input_fastqs, bwa_index, bam_sort, ch_fasta )
      
    }
}


