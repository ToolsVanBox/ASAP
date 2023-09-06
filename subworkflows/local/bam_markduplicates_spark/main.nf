//
// BAM MARKDUPLICATES

 include { GATK4_MARKDUPLICATES_SPARK } from '../../../modules/nf-core/gatk4/markduplicatesspark/main'        

include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_SPARK {
  take:
    ch_bams   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ path(fai) ]
    ch_dict // channel: [ path(dict) ]

  main:
    ch_versions = Channel.empty()
    
    ch_fasta_2 = ch_fasta.map{ meta, fasta -> [ fasta ] }

    GATK4_MARKDUPLICATES_SPARK( ch_bams, ch_fasta_2, ch_fai, ch_dict )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.versions)

    SAMTOOLS_INDEX ( GATK4_MARKDUPLICATES_SPARK.out.output )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = GATK4_MARKDUPLICATES_SPARK.out.output
      .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
      .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
      .map {
          meta, bam, bai, csi ->
              if (bai) {
                  [ meta, bam, bai ]
              } else {
                  [ meta, bam, csi ]
              }
      }
    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
  
  emit:
    bam      = GATK4_MARKDUPLICATES_SPARK.out.output     // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), path(csi) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}


