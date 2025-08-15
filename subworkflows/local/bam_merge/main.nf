//
// Merge bams files with Samtools
//

// Include nf-core modules
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'

// Include nf-core subworkflows
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow BAM_MERGE {
    take:
    ch_bams        // channel (mandatory): [ val(meta), [ path(bams) ] ]
    ch_fasta       // channel (optional) : [ val(meta2), path(fasta) ]
    ch_fai         // channel (optional) : [ val(meta3), path(fai) ]
    
    main:
        ch_versions = Channel.empty()
        ch_merged_bam = Channel.empty()
        ch_merged_bai = Channel.empty()
        ch_merged_csi = Channel.empty()
        ch_merged_stats = Channel.empty()
        ch_merged_flagstat = Channel.empty()
        ch_merged_idxstats = Channel.empty()
            
        for ( tool in params.bam_merge.tool ) {
            def tool = tool.toLowerCase()      
            def known_tool = false    

            if ( tool == "samtools" ) {
                def ch_bams_samtools = ch_bams
                     .map{ meta, bam ->              
                        meta = meta + [ merge: "samtools" ]
                        meta.id = meta.id+".samtools"
                        [ meta, bam ]
                    }
                // Run Samtools merge
                SAMTOOLS_MERGE ( ch_bams_samtools, ch_fasta, ch_fai )
                ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

                // Sort, index BAM file and run samtools stats, flagstat and idxstats    
                BAM_SORT_STATS_SAMTOOLS ( SAMTOOLS_MERGE.out.bam, ch_fasta )
                ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

                ch_merged_bam = ch_merged_bam.mix( BAM_SORT_STATS_SAMTOOLS.out.bam )
                ch_merged_bai = ch_merged_bai.mix( BAM_SORT_STATS_SAMTOOLS.out.bai )
                ch_merged_csi = ch_merged_csi.mix( BAM_SORT_STATS_SAMTOOLS.out.csi )
                ch_merged_stats = ch_merged_stats.mix( BAM_SORT_STATS_SAMTOOLS.out.stats )
                ch_merged_flagstat = ch_merged_flagstat.mix( BAM_SORT_STATS_SAMTOOLS.out.flagstat )
                ch_merged_idxstats = ch_merged_idxstats.mix( BAM_SORT_STATS_SAMTOOLS.out.idxstats )

                known_tool = true
            }

            if ( ! known_tool ) {
                println ("WARNING: Skip ${tool}, because it's not known as a merge bam tool, this tool is not build in (yet).")
            }
        }

    emit:
        bam      = ch_merged_bam      // channel: [ val(meta), path(bam) ]
        bai      = ch_merged_bai      // channel: [ val(meta), path(bai) ]
        csi      = ch_merged_csi      // channel: [ val(meta), path(csi) ]
        stats    = ch_merged_stats    // channel: [ val(meta), path(stats) ]
        flagstat = ch_merged_flagstat // channel: [ val(meta), path(flagstat) ]
        idxstats = ch_merged_idxstats // channel: [ val(meta), path(idxstats) ]

        versions = ch_versions                          // channel: [ path(versions.yml) ]
}