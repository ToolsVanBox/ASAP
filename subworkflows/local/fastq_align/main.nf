//
// FASTQC

include { FASTQ_ALIGN_BWA } from '../../nf-core/fastq_align_bwa/main.nf'
include { PREPARE_GENOME } from '../prepare_genome/main.nf'
include { BAM_MERGE } from '../bam_merge/main.nf'

workflow FASTQ_ALIGN {
  take:
    input_fastqs
  main:

      // Convert input tools to UpperCase
    fastqc_align_tools = params.fastq_align.tool.toString().toUpperCase()
    genome_fasta = params.genomes[params.genome].fasta

    if ( fastqc_align_tools == "BWA_MEM" | fastqc_align_tools.contains("BWA_MEM") ) {
        tool = "BWA_MEM"
        
        ch_fasta = Channel.value( file(genome_fasta) )
          .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }
        
        bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
        
        if ( !bwa_index_folder.exists() ) {
          PREPARE_GENOME( ch_fasta, tool )

          bwa_index = PREPARE_GENOME.out.bwa
            .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
          
        } else {
          bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
            .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
        }

        bam_sort = params.fastq_align.bam_sort
        FASTQ_ALIGN_BWA( input_fastqs, bwa_index, bam_sort, ch_fasta )
        
        bwa_bams = FASTQ_ALIGN_BWA.out.bam
          .map{ meta, bam -> [ [id:meta.sample_id, tool:"${tool}"], bam ] }
          .groupTuple()
    } else{
        error( "Cannot align fastq files with tool ${fastqc_align_tools}. This module is not build in (yet)." )
    }

    // Merge bam files 
    input_bam_merge = Channel.empty()
    input_bam_merge = input_bam_merge.mix( bwa_bams )
    
    input_bam_merge.view()
    ch_fai = Channel.value( file(genome_fasta+".fai") )
      .map{ genome_fasta -> [ [ id:'fai' ], genome_fasta ] }

    BAM_MERGE( input_bam_merge, ch_fasta, ch_fai )    
}


