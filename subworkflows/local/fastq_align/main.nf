//
// FASTQC

include { FASTQ_ALIGN_BWA } from '../../nf-core/fastq_align_bwa/main.nf'
include { PREPARE_GENOME } from '../prepare_genome/main.nf'
include { BAM_MERGE } from '../bam_merge/main.nf'

include { FASTQ_ALIGN_BWAMEM2 } from '../fastq_align_bwamem2/main.nf'

workflow FASTQ_ALIGN {
  take:
    input_fastqs
  main:

      // Convert input tools to UpperCase
    fastq_align_tools = params.fastq_align.tool.toString().toUpperCase()

    genome_fasta = params.genomes[params.genome].fasta
    ch_fasta = Channel.value( file(genome_fasta) )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( file(genome_fasta+".fai") )
      .map{ genome_fasta -> [ [ id:'fai' ], genome_fasta ] }

    input_bam_merge = Channel.empty()    

    if ( fastq_align_tools == "BWA_MEM1" | fastq_align_tools.contains("BWA_MEM1") ) {
        tool = "bwa"
                
        bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
        
        if ( !bwa_index_folder.exists() ) {
          PREPARE_GENOME( ch_fasta, tool )

          bwa_index = PREPARE_GENOME.out.index
            .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
          
        } else {
          bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
            .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
        }

        bam_sort = params.fastq_align.bam_sort
        FASTQ_ALIGN_BWA( input_fastqs, bwa_index, bam_sort, ch_fasta )
        
        bwamem1_bams = FASTQ_ALIGN_BWA.out.bam
          .map{ meta, bam -> [ [sample_id: meta.sample_id, id:meta.sample_id+".${tool}", tool:"${tool}"], bam ] }
          .groupTuple()

        input_bam_merge = input_bam_merge.mix( bwamem1_bams )
    
    } 

    if ( fastq_align_tools == "BWA_MEM2" | fastq_align_tools.contains("BWA_MEM2") ) {
        tool = "bwamem2"
                
        bwamem2_index_folder = new File("${params.genomes[params.genome].bwamem2}")
        
        if ( !bwamem2_index_folder.exists() ) {
          PREPARE_GENOME( ch_fasta, tool )

          bwamem2_index = PREPARE_GENOME.out.index
            .map{ bwamem2_index -> [ [ id:'index' ], bwamem2_index ] }
          
        } else {
          bwamem2_index = Channel.value( file("${params.genomes[params.genome].bwamem2}/*") )
            .map{ bwamem2_index -> [ [ id:'index' ], bwamem2_index ] }
        }

        bam_sort = params.fastq_align.bam_sort
        FASTQ_ALIGN_BWAMEM2( input_fastqs, bwamem2_index, bam_sort, ch_fasta )
        
        bwamem2_bams = FASTQ_ALIGN_BWAMEM2.out.bam
          .map{ meta, bam -> [ [sample_id: meta.sample_id, id:meta.sample_id+".${tool}", tool:"${tool}"], bam ] }
          .groupTuple()

        input_bam_merge = input_bam_merge.mix( bwamem2_bams )
    }

    BAM_MERGE( input_bam_merge, ch_fasta, ch_fai )

    bams = BAM_MERGE.out.bam
      .join( BAM_MERGE.out.bai )

    emit:
      bam = bams
}


