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
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ] 
    ch_fasta // channel : [ meta, path(fasta) ]
    ch_fai // channel: [ meta, path(fai) ]
    ch_dict // channel: [ meta, path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_dedup_bam = Channel.empty()
    ch_dedup_bai = Channel.empty()
      
    for ( tool in params.bam_markduplicates.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      // Run Picard Markduplicates
      if ( tool == "picard" ) {          
          ch_bams_picard = ch_bam_bai
            .map{ meta, bam, bai -> 
              meta = meta + [ dedup: "picard" ]
              meta.id = meta.id+".picard.dedup"
              [ meta, bam ]
            }
          
          BAM_MARKDUPLICATES_PICARD( ch_bams_picard, ch_fasta, ch_fai )
          ch_versions = ch_versions.mix( BAM_MARKDUPLICATES_PICARD.out.versions )
          
          ch_dedup_bam = ch_dedup_bam.mix( BAM_MARKDUPLICATES_PICARD.out.bam )
          ch_dedup_bai = ch_dedup_bai.mix( BAM_MARKDUPLICATES_PICARD.out.bai )

          known_tool = true
      } 

      // Run Samtools Markduplicates
      if ( tool == "samtools" ) {   
          ch_bams_samtools = ch_bam_bai
            .map{ meta, bam, bai -> 
              meta = meta + [ dedup: "samtools" ]
              meta.id = meta.id+".samtools.dedup"
              [ meta, bam ]
            }
        
          BAM_MARKDUPLICATES_SAMTOOLS( ch_bams_samtools, ch_fasta )
          ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions)

          ch_dedup_bam = ch_dedup_bam.mix( BAM_MARKDUPLICATES_SAMTOOLS.out.bam )
          ch_dedup_bai = ch_dedup_bai.mix( BAM_MARKDUPLICATES_SAMTOOLS.out.bai )
          
          known_tool = true

      }

      // Run GATK4 Spark Markduplicates
      if ( tool == "gatk4spark" ) {
          ch_bams_spark = ch_bam_bai
            .map{ meta, bam, bai -> 
              meta = meta + [ dedup: "gatk4spark" ]
              meta.id = meta.id+".gatk4spark.dedup"
              [ meta, bam ]
            }

          BAM_MARKDUPLICATES_GATK4SPARK( ch_bams_spark, ch_fasta, ch_fai, ch_dict )
          ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_GATK4SPARK.out.versions)

          ch_dedup_bam = ch_dedup_bam.mix( BAM_MARKDUPLICATES_GATK4SPARK.out.bam )
          ch_dedup_bai = ch_dedup_bai.mix( BAM_MARKDUPLICATES_GATK4SPARK.out.bai )

          known_tool = true

        }

      // Run GATK4 Markduplicates
      if ( tool == "gatk4" ) {
        ch_bams_gatk = ch_bam_bai
            .map{ meta, bam, bai -> 
              meta = meta + [ dedup: "gatk4" ]
              meta.id = meta.id+".gatk4.dedup"
              [ meta, bam ]
            }

        BAM_MARKDUPLICATES_GATK4( ch_bams_gatk, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_MARKDUPLICATES_GATK4.out.versions )

        ch_dedup_bam = ch_dedup_bam.mix( BAM_MARKDUPLICATES_GATK4.out.bam )
        ch_dedup_bai = ch_dedup_bai.mix( BAM_MARKDUPLICATES_GATK4.out.bai )

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam markdup tool, this tool is not build in (yet).")
      }
    }
    
  emit:
    bam = ch_dedup_bam // channel: [ meta, bam ]
    bai = ch_dedup_bai // channel: [ meta, bai ]
    versions = ch_versions // channel: [ versions.yml ]

}