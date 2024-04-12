//
// BAM MARKDUPLICATES SAMTOOLS

// Include nf-core subworkflows
include { BAM_MARKDUPLICATES_SAMTOOLS as BAM_MARKDUPLICATES_SAMTOOLS_NFCORE } from '../../../subworkflows/nf-core/bam_markduplicates_samtools/main'                         
include { BAM_STATS_SAMTOOLS    } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

// Include nf-core modules
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MARKDUPLICATES_SAMTOOLS {
  take:
    ch_bams   // channel: [ val(meta), path(bam) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    
  main:
    ch_versions = Channel.empty()
    
    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    
    BAM_MARKDUPLICATES_SAMTOOLS_NFCORE( ch_bams, fasta )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS_NFCORE.out.versions)

    SAMTOOLS_INDEX ( BAM_MARKDUPLICATES_SAMTOOLS_NFCORE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = BAM_MARKDUPLICATES_SAMTOOLS_NFCORE.out.bam
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
    bam      = BAM_MARKDUPLICATES_SAMTOOLS_NFCORE.out.bam     // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai            // channel: [ val(meta), path(bai) ]
    csi      = SAMTOOLS_INDEX.out.csi            // channel: [ val(meta), path(csi) ]

    stats    = BAM_STATS_SAMTOOLS.out.stats      // channel: [ val(meta), path(stats) ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat   // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats   // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                       // channel: [ versions.yml ]
}


