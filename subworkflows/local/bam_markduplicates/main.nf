//
// BAM MARKDUPLICATES

include { BAM_MARKDUPLICATES_PICARD } from '../../nf-core/bam_markduplicates_picard/main'                             
include { BAM_MARKDUPLICATES_SAMTOOLS } from '../../local/bam_markduplicates_samtools/main'        
include { BAM_MARKDUPLICATES_SPARK } from '../../local/bam_markduplicates_spark/main'       
include { BAM_MARKDUPLICATES_GATK } from '../../local/bam_markduplicates_gatk/main'


workflow BAM_MARKDUPLICATES {
  take:
    input_bams 
  main:
    
    // Convert input tools to UpperCase
    bam_markduplicates_tools = params.bam_markduplicates.tool.toString().toUpperCase()

    ch_bams = input_bams.map{ meta, bam, bai -> [ meta + [ step:"bam_markduplicates" ], bam ] } 
    
    def fasta_file = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai_file = file( fasta_file.toString()+".fai", checkIfExists: true )
    def dict_file = file( fasta_file.toString()+".dict", checkIfExists: true )

    ch_fasta = Channel.value( fasta_file )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( fai_file )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }
    ch_dict = Channel.value( dict_file )

    ch_fasta_2 = ch_fasta.map{ meta, fasta -> [ fasta ] }
    ch_fai_2 = ch_fai.map{ meta, fai -> [ fai ] }
    
    ch_versions = Channel.empty()

    bams = Channel.empty()

    if ( bam_markduplicates_tools == "PICARD" | bam_markduplicates_tools.contains("PICARD") ) {
        tool = "picard"

        ch_bams = ch_bams
          .map{ meta, bam -> 
            meta.id = meta.id+".${tool}.dedup" 
            [ meta, bam ]
        }

        BAM_MARKDUPLICATES_PICARD( ch_bams, ch_fasta, ch_fai )   

        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
        
        bams = bams.mix(BAM_MARKDUPLICATES_PICARD.out.bam.join( BAM_MARKDUPLICATES_PICARD.out.bai) )

    } 

    if ( bam_markduplicates_tools == "SAMTOOLS" | bam_markduplicates_tools.contains("SAMTOOLS") ) {
        tool = "samtools"          

        ch_bams = ch_bams
          .map{ meta, bam -> 
            meta.id = meta.id+".${tool}.dedup" 
            [ meta, bam ]
        }
        
        BAM_MARKDUPLICATES_SAMTOOLS( ch_bams, ch_fasta_2 )
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions)

        bams = bams.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join( BAM_MARKDUPLICATES_SAMTOOLS.out.bai) )

    }

    if ( bam_markduplicates_tools == "SPARK" | bam_markduplicates_tools.contains("SPARK") ) {
      tool = "spark"

      ch_bams = ch_bams
        .map{ meta, bam -> 
          meta.id = meta.id+".${tool}.dedup" 
          [ meta, bam ]
      }      
      BAM_MARKDUPLICATES_SPARK( ch_bams, ch_fasta, ch_fai_2, ch_dict )
      ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SPARK.out.versions)

      bams = bams.mix(BAM_MARKDUPLICATES_SPARK.out.bam.join( BAM_MARKDUPLICATES_SPARK.out.bai) )

    }

    if ( bam_markduplicates_tools == "GATK" | bam_markduplicates_tools.contains("GATK") ) {
      tool = "gatk"

      ch_bams = ch_bams
        .map{ meta, bam -> 
          meta.id = meta.id+".${tool}.dedup" 
          [ meta, bam ]
      }      
      BAM_MARKDUPLICATES_GATK( ch_bams, ch_fasta, ch_fai_2 )
      ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_GATK.out.versions)

      bams = bams.mix(BAM_MARKDUPLICATES_GATK.out.bam.join( BAM_MARKDUPLICATES_GATK.out.bai) )

    }

  emit:
    bam = bams

}


