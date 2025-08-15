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
    ch_fasta // channel: [ meta, path(fasta) ]
  main:
    ch_versions = Channel.empty()
    ch_align_bam = Channel.empty()
    ch_align_bai = Channel.empty()    

    for ( tool in params.fastq_align.tool ) {
      tool = tool.toLowerCase()      
      def known_tool = false    

      // Run BWA MEM 1
      if ( tool == "bwamem1" ) {
         
          def bwa_index_folder = file( params.genomes[params.genome].bwa, checkIfExists: false )
          
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

          ch_fastq_bwa = ch_fastq
            .map{ meta, bam ->              
              meta = meta + [align: "bwamem1" ]
              meta.id = meta.id+".bwamem1"
              [ meta, bam ]
            }
          
          FASTQ_ALIGN_BWA( ch_fastq_bwa, bwa_index, bam_sort, ch_fasta )
          ch_versions = ch_versions.mix( FASTQ_ALIGN_BWA.out.versions.first() )
          
          ch_align_bam = ch_align_bam.mix( FASTQ_ALIGN_BWA.out.bam )
          ch_align_bai = ch_align_bam.mix( FASTQ_ALIGN_BWA.out.bai )

          known_tool = true
      } 

      // Run BWA MEM 2
      if ( tool == "bwamem2" ) {
                    
          def bwamem2_index_folder = file( params.genomes[params.genome].bwamem2, checkIfExists: false )

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

          ch_fastq_bwamem2 = ch_fastq
            .map{ meta, bam ->
              meta = meta + [ align: "bwamem2" ]
              meta.id = meta.id+".bwamem2"
              [ meta, bam ]
            }

          FASTQ_ALIGN_BWAMEM2( ch_fastq_bwamem2, bwamem2_index, bam_sort, ch_fasta )
          ch_versions = ch_versions.mix( FASTQ_ALIGN_BWAMEM2.out.versions.first() )

          ch_align_bam = ch_align_bam.mix( FASTQ_ALIGN_BWAMEM2.out.bam )
          ch_align_bai = ch_align_bam.mix( FASTQ_ALIGN_BWAMEM2.out.bai )
          
          known_tool = true
      }

      if ( ! known_tool ) {
        println ("WARNING: Skip ${tool}, because it's not known as a Fastq alinger tool, this tool is not build in (yet).")
      }
    }

    emit:
      bam = ch_align_bam // channel: [ meta, bam ]
      bai = ch_align_bai // channel: [ meta, bai ]
      versions = ch_versions // channel: [ versions.yml ]

}