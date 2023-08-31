//
// Merge bams files with Samtools
//

include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'                                                                                       
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow BAM_MERGE {
    take:
    ch_bams        // channel (mandatory): [ val(meta), [ path(bams) ] ]
    ch_fasta       // channel (optional) : [ val(meta2), path(fasta) ]
    ch_fai         // channel (optional) : [ val(meta3), path(fai) ]
    
    main:
    ch_versions = Channel.empty()

    
    // Convert input tools to UpperCase
    bam_merge_tools = params.bam_merge.tool.toString().toUpperCase()
    
    if ( bam_merge_tools == "SAMTOOLS" | bam_merge_tools.contains("SAMTOOLS") ) {
        // Merge bam file with samtools
        SAMTOOLS_MERGE ( ch_bams, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())
    } else {
        error( "Cannot merge bam file with tool ${bam_merge_tools}. This module is not build in (yet)." )
    }
    // Sort, index BAM file and run samtools stats, flagstat and idxstats    

    BAM_SORT_STATS_SAMTOOLS ( SAMTOOLS_MERGE.out.bam, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
