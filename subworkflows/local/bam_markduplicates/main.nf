//
// BAM MARKDUPLICATES

// Include nf-core subworkflows
include { BAM_MARKDUPLICATES_PICARD } from '../../../subworkflows/nf-core/bam_markduplicates_picard/main'                             

// Include local subworkflows
include { BAM_MARKDUPLICATES_SAMTOOLS } from '../../../subworkflows/local/bam_markduplicates_samtools/main'        
include { BAM_MARKDUPLICATES_GATK4SPARK } from '../../../subworkflows/local/bam_markduplicates_gatk4spark/main'       
include { BAM_MARKDUPLICATES_GATK4 } from '../../../subworkflows/local/bam_markduplicates_gatk4/main'


workflow BAM_MARKDUPLICATES {
  take:
    ch_bams // channel: [ meta, path(bam) ] 
  main:
    ch_versions = Channel.empty()
    ch_dedup_bams = Channel.empty()

    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( fasta.toString()+".fai", checkIfExists: true )
    def fasta_dict = file( fasta.toString()+".dict", checkIfExists: true )
    def dict = file( fasta.toString().replace(".fasta",".dict" ), checkIfExists: true )
    
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( fai )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }

    ch_bams = ch_bams.map{ meta, bam, bai -> [ meta, bam ] }
        
    for ( tool in params.bam_markduplicates.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Markduplicates
      if ( tool == "picard" ) {          
          ch_bams_picard = ch_bams.map{ meta, bam -> [ [ sample_id: meta.sample_id, id: meta.id+".picard.dedup" ], bam ] }
          
          BAM_MARKDUPLICATES_PICARD( ch_bams_picard, ch_fasta, ch_fai )
          ch_versions = ch_versions.mix( BAM_MARKDUPLICATES_PICARD.out.versions )
          
          ch_dedup_bams = ch_dedup_bams.mix( BAM_MARKDUPLICATES_PICARD.out.bam.join( BAM_MARKDUPLICATES_PICARD.out.bai ) )

          known_tool = true
      } 

      // Run Samtools Markduplicates
      if ( tool == "samtools" ) {   
          ch_bams_samtools = ch_bams.map{ meta, bam -> [ [ sample_id: meta.sample_id, id: meta.id+".samtools.dedup" ], bam ] }
        
          BAM_MARKDUPLICATES_SAMTOOLS( ch_bams_samtools, fasta )
          ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions)

          ch_dedup_bams = ch_dedup_bams.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join( BAM_MARKDUPLICATES_SAMTOOLS.out.bai) )
          
          known_tool = true

      }

      // Run GATK4 Spark Markduplicates
      if ( tool == "gatk4spark" ) {
          ch_bams_spark = ch_bams.map{ meta, bam -> [ [ sample_id: meta.sample_id, id: meta.id+".gatk4spark.dedup" ], bam ] }

          BAM_MARKDUPLICATES_GATK4SPARK( ch_bams_spark, fasta, fai, dict )
          ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_GATK4SPARK.out.versions)

          ch_dedup_bams = ch_dedup_bams.mix(BAM_MARKDUPLICATES_GATK4SPARK.out.bam.join( BAM_MARKDUPLICATES_GATK4SPARK.out.bai) )

          known_tool = true

        }

      // Run GATK4 Markduplicates
      if ( tool == "gatk4" ) {        
        ch_bams_gatk = ch_bams.map{ meta, bam -> [ [ sample_id: meta.sample_id, id: meta.id+".gatk4.dedup" ], bam ] } 

        BAM_MARKDUPLICATES_GATK4( ch_bams_gatk, fasta, fai )
        ch_versions = ch_versions.mix( BAM_MARKDUPLICATES_GATK4.out.versions )

        ch_dedup_bams = ch_dedup_bams.mix(BAM_MARKDUPLICATES_GATK4.out.bam.join( BAM_MARKDUPLICATES_GATK4.out.bai ) )

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam markdup tool, this tool is not build in (yet).")
      }
    }
  emit:
    bam = ch_dedup_bams // channel: [ meta, bam, bai ]
    versions = ch_versions // channel: [ versions.yml ]

}


