//
// FASTQC

// Include nf-core subworkflows
include { FASTQ_ALIGN_BWA } from '../../../subworkflows/nf-core/fastq_align_bwa/main.nf'

// Include local subworkflows
include { PREPARE_GENOME } from '../../../subworkflows/local/prepare_genome/main.nf'
include { BAM_MERGE } from '../../../subworkflows/local/bam_merge/main.nf'
include { FASTQ_ALIGN_BWAMEM2 } from '../../../subworkflows/local/fastq_align_bwamem2/main.nf'

workflow FASTQ_ALIGN {
  take:
    ch_fastq // channel: [ meta, path(fastq) ]
  main:
    ch_versions = Channel.empty()
    ch_bams_to_merge = Channel.empty()    

    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( fasta.toString()+".fai", checkIfExists: true )
    
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( fai )
      .map{ genome_fasta -> [ [ id:'fai' ], genome_fasta ] }

        
    for ( tool in params.fastq_align.tool ) {
      tool = tool.toLowerCase()      
      def known_tool = false    

      // Run BWA MEM 1
      if ( tool == "bwamem1" ) {
         
          bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
          
          // Create BWA index if it does not exists
          if ( !bwa_index_folder.exists() ) {
            PREPARE_GENOME( ch_fasta, 'bwamem1' )
            ch_versions = ch_versions.mix( PREPARE_GENOME.out.versions.first() )

            bwa_index = PREPARE_GENOME.out.index
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
            
          } else {
            bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }
          }

          // Align fastq per lane
          bam_sort = params.fastq_align.bam_sort
          FASTQ_ALIGN_BWA( ch_fastq, bwa_index, bam_sort, ch_fasta )
          ch_versions = ch_versions.mix( FASTQ_ALIGN_BWA.out.versions.first() )
          
          ch_bwamem1_bams = FASTQ_ALIGN_BWA.out.bam
            .map{ meta, bam ->
              [ [sample_id: meta.sample_id, id:meta.sample_id+".bwamem1"], bam ] 
            }
            .groupTuple()

          ch_bams_to_merge = ch_bams_to_merge.mix( ch_bwamem1_bams )

          known_tool = true
      } 

      // Run BWA MEM 2
      if ( tool == "bwamem2" ) {
                    
          bwamem2_index_folder = new File("${params.genomes[params.genome].bwamem2}")
          
          // Create BWA MEM 2 index if it does not exists
          if ( !bwamem2_index_folder.exists() ) {
            PREPARE_GENOME( ch_fasta, 'bwamem2' )
            ch_versions = ch_versions.mix( PREPARE_GENOME.out.versions.first() )

            bwamem2_index = PREPARE_GENOME.out.index
              .map{ bwamem2_index -> [ [ id:'index' ], bwamem2_index ] }
            
          } else {
            bwamem2_index = Channel.value( file("${params.genomes[params.genome].bwamem2}/*") )
              .map{ bwamem2_index -> [ [ id:'index' ], bwamem2_index ] }
          }

          // Align fastq per lane
          bam_sort = params.fastq_align.bam_sort
          FASTQ_ALIGN_BWAMEM2( ch_fastq, bwamem2_index, bam_sort, ch_fasta )
          ch_versions = ch_versions.mix( FASTQ_ALIGN_BWAMEM2.out.versions.first() )

          ch_bwamem2_bams = FASTQ_ALIGN_BWAMEM2.out.bam
            .map{ meta, bam -> 
              [ [sample_id: meta.sample_id, id:meta.sample_id+".bwamem2"], bam ] 
            }
            .groupTuple()

          ch_bams_to_merge = ch_bams_to_merge.mix( ch_bwamem2_bams )
          
          known_tool = true
      }

      if ( ! known_tool ) {
        println ("WARNING: Skip ${tool}, because it's not known as a Fastq alinger tool, this tool is not build in (yet).")
      }
    }

    // Merge bam files
    BAM_MERGE( ch_bams_to_merge, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix( BAM_MERGE.out.versions.first() )

    ch_bams = BAM_MERGE.out.bam
      .join( BAM_MERGE.out.bai )

    emit:
      bam = ch_bams // channel: [ meta, bam, bai ]
      versions = ch_versions // channel: [ versions.yml ]

}


