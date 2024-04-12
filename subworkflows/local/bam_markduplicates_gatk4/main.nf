//
// BAM MARKDUPLICATES GATK4

// Include nf-core modules
include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'        
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'

// Include nf-core subworkflows
include { BAM_STATS_SAMTOOLS    } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_GATK4 {
  take:
    ch_bams   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(metae), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

  main:
    ch_versions = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    fai = ch_fai.map{ meta, fai -> [ fai ] }

    GATK4_MARKDUPLICATES( ch_bams, fasta, fai )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    SAMTOOLS_INDEX ( GATK4_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = GATK4_MARKDUPLICATES.out.bam
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
    bam      = GATK4_MARKDUPLICATES.out.bam     // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), path(csi) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}


