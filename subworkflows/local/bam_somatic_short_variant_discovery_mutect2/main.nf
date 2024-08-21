//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MUTECT2                   as MUTECT2 }                   from '../../../modules/local/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR }  from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL}  from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL} from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR } from '../../../modules/nf-core/gatk4/gatherpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION }    from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS }         from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_MERGEVCFS                 as MERGE_MUTECT2               } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS            } from '../../../modules/nf-core/gatk4/mergemutectstats/main'



workflow BAM_SOMATIC_SHORT_VARIANT_DISCOVERY_MUTECT2 {
    take:
        ch_bam_bai_tumor // channel: [ val(meta), path(bam), path(bai) ]
        ch_bam_bai_normal // channel: [ val(meta), path(bam), path(bai) ]
        ch_split_intervals // channel: [ intervals ]
        ch_fasta // channel: [ meta, path(fasta) ]
        ch_fai // channel: [ meta, path(fai) ]
        ch_dict // channel: [ meta, path(dict) ]

    main:
        ch_versions = Channel.empty()
        
        def germline_resource = file( params.genomes[params.genome].somatic_germline_resource, checkIfExists: true )
        def germline_resource_tbi = file( germline_resource.toString()+".tbi", checkIfExists: true )
        def pon = file( params.genomes[params.genome].somatic_pon, checkIfExists: true )
        def pon_tbi = file( pon.toString()+".tbi", checkIfExists: true )
        def dict = ch_dict.map{ meta, dict -> [ dict ] }
        
        ch_normal_ids = ch_bam_bai_normal
            .map{ meta, bam, bai ->
            [ [id: meta.run_id ], meta.sample ]
            }
            .groupTuple()

        // Group tumor and normal samples while adding meta info 
        ch_input_mutect2 = ch_bam_bai_tumor
            .map{ meta, bam, bai ->
            [ [ id: meta.run_id, run_id: meta.run_id ], bam, bai ]
            }
            .mix( ch_bam_bai_normal
                .map{ meta, bam, bai ->             
                [ [ id: meta.run_id, run_id: meta.run_id ], bam, bai ]
                }
            )
            .groupTuple()
            .combine( ch_split_intervals )
            .map { meta, bam, bai, meta2, interval_file ->
                meta = meta + [ interval: meta2.id ]
                meta.id = meta.id+".mutect2."+meta2.id
                [ meta, bam, bai, interval_file ]
            }
            .combine( ch_normal_ids )
            .map{ meta, bam, bai, interval, meta2, normal_ids ->
                meta = meta + [ normal_ids: normal_ids ]
                meta = meta + [ calling_type: "somatic" ] 
                meta = meta + [ short_variant_caller: "mutect2" ]
                [ meta, bam, bai, interval ]
            }

        // Run mutect2 and merge output 
        MUTECT2( ch_input_mutect2, ch_fasta, ch_fai, ch_dict, germline_resource, germline_resource_tbi, pon, pon_tbi )
        ch_versions = ch_versions.mix(MUTECT2.out.versions)
        
        ch_merge_mutect2 = MUTECT2.out.vcf
            .map{ meta, vcf ->
                meta = meta - meta.subMap("id","interval","normal_ids")
                meta = meta + [ id: meta.run_id ]
                [ meta, vcf ]
            }
            .groupTuple()

        MERGE_MUTECT2(ch_merge_mutect2, ch_dict)
        ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)

        ch_merge_mutect2_stats = MUTECT2.out.stats
            .map{ meta, stats ->
                meta = meta - meta.subMap("id","interval","normal_ids")
                meta = meta + [ id: meta.run_id ]
                [ meta, stats ]
            }
            .groupTuple()

        MERGEMUTECTSTATS(ch_merge_mutect2_stats)   
        ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)

        // Get all the pileups 
        ch_pileup_tumor_input = ch_bam_bai_tumor
            .combine( ch_split_intervals )
            .map { meta, bam, bai, meta2, interval_file ->
              meta = meta + [ tumor_id: meta.id ]
              meta = meta + [ interval: meta2.id ]
              meta.id = meta.id+".pileupsummary."+meta2.id
              [ meta, bam, bai, interval_file ]
            }

        ch_pileup_normal_input = ch_bam_bai_normal
            .combine( ch_split_intervals )
            .map { meta, bam, bai, meta2, interval_file ->
              meta = meta + [ normal_id: meta.id ]
              meta = meta + [ interval: meta2.id ]
              meta.id = meta.id+".pileupsummary."+meta2.id
              [ meta, bam, bai, interval_file ]
            }

        GETPILEUPSUMMARIES_TUMOR( ch_pileup_tumor_input, ch_fasta, ch_fai, ch_dict, germline_resource, germline_resource_tbi )
        GETPILEUPSUMMARIES_NORMAL( ch_pileup_normal_input, ch_fasta, ch_fai, ch_dict, germline_resource, germline_resource_tbi )
        ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions.first())
        ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions.first())

        ch_gather_pileup_sum_tumor = GETPILEUPSUMMARIES_TUMOR.out.table
            .map{ meta, table ->
              meta = meta - meta.subMap("interval")
              meta.id = meta.id.replaceAll(/.\d+$/, "")
              [ meta, table ]
            }
            .groupTuple()

        ch_gather_pileup_sum_normal = GETPILEUPSUMMARIES_NORMAL.out.table
            .map{ meta, table ->
              meta = meta - meta.subMap("interval")
              meta.id = meta.id.replaceAll(/.\d+$/, "")
              [ meta, table ]
            }
            .groupTuple()

        GATHERPILEUPSUMMARIES_TUMOR( ch_gather_pileup_sum_tumor, dict ) 
        GATHERPILEUPSUMMARIES_NORMAL( ch_gather_pileup_sum_normal, dict )  
        ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_TUMOR.out.versions)
        ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
        
        // Learn the orientation model
        ch_learnorientationmodel = MUTECT2.out.f1r2
            .map{ meta, f1r2_file -> 
                meta = meta - meta.subMap("id","interval","normal_ids")
                meta = meta + [ id: meta.run_id ]
                [ meta, f1r2_file ]
            }
            .groupTuple()
            
        LEARNREADORIENTATIONMODEL( ch_learnorientationmodel )
        ch_versions = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)

        ch_learnorientationmodel_artifactprior = LEARNREADORIENTATIONMODEL.out.artifactprior
            .map{ meta, arti -> 
                meta = meta - meta.subMap("id")
                meta = meta + [ id: meta.run_id ]
                [meta, arti]
            }

        ch_pileup_tumor = GATHERPILEUPSUMMARIES_TUMOR.out.table
        ch_pileup_normal = GATHERPILEUPSUMMARIES_NORMAL.out.table

        // Calculate the contamination in the sample
        ch_calculate_cont = ch_pileup_tumor
            .combine( ch_pileup_normal )
            .map{ meta, tumor_table, meta2, normal_table ->
                meta = meta - meta.subMap("id","sample_type")
                meta = meta + [ id: meta.sample+"-"+meta2.sample ]                
                meta = meta + [ normal_id: meta2.normal_id ]
                meta = meta - meta.subMap("sample")
                [ meta, tumor_table, normal_table ]
            }

        CALCULATECONTAMINATION ( ch_calculate_cont )
        ch_versions   = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

        ch_segmentation_output = CALCULATECONTAMINATION.out.segmentation
            .map{ meta, segm ->
                [ [id: meta.run_id, run_id: meta.run_id, calling_type: "somatic", short_variant_caller: "mutect2" ], segm ]
            }
            .groupTuple()
        
        ch_contamination_output = CALCULATECONTAMINATION.out.contamination
            .map{ meta, cont ->
                [ [id: meta.run_id, run_id: meta.run_id, calling_type: "somatic", short_variant_caller: "mutect2" ], cont ]
            }
            .groupTuple()
        
        // Filter mutect vcf file 
        ch_filter_mutect = MERGE_MUTECT2.out.vcf
            .join( MERGE_MUTECT2.out.tbi )
            .join( MERGEMUTECTSTATS.out.stats )
            .join( ch_learnorientationmodel_artifactprior )
            .join( ch_segmentation_output )
            .join( ch_contamination_output )

        ch_filter_mutect2 = ch_filter_mutect.map{ meta, vcf, tbi, stats, orientation, seg, cont -> 
        meta.id = meta.id+".filt_mutect2"
        [ meta, vcf, tbi, stats, orientation, seg, cont, [] ] }

        FILTERMUTECTCALLS(ch_filter_mutect2, ch_fasta, ch_fai, ch_dict)
        ch_versions   = ch_versions.mix(FILTERMUTECTCALLS.out.versions)

    emit:
        mutect2_vcf            = MUTECT2.out.vcf.collect()                             // channel: [ val(meta), path(vcf) ]
        mutect2_tbi            = MUTECT2.out.tbi.collect()                             // channel: [ val(meta), path(tbi) ]
        mutect2_stats          = MUTECT2.out.stats.collect()                           // channel: [ val(meta), path(stats) ]
        mutect2_f1r2           = MUTECT2.out.f1r2.collect()                            // channel: [ val(meta), path(f1r2) ]

        artifact_priors        = LEARNREADORIENTATIONMODEL.out.artifactprior.collect() // channel: [ val(meta), path(artifactprior) ]

        pileup_table_tumor     = GETPILEUPSUMMARIES_TUMOR.out.table.collect()          // channel: [ val(meta), path(table) ]
        pileup_table_normal    = GETPILEUPSUMMARIES_NORMAL.out.table.collect()         // channel: [ val(meta), path(table) ]

        contamination_table    = CALCULATECONTAMINATION.out.contamination.collect()    // channel: [ val(meta), path(table) ]
        segmentation_table     = CALCULATECONTAMINATION.out.segmentation.collect()     // channel: [ val(meta), path(table) ]

        filtered_vcf           = FILTERMUTECTCALLS.out.vcf.collect()                   // channel: [ val(meta), path(vcf) ]
        filtered_tbi           = FILTERMUTECTCALLS.out.tbi.collect()                   // channel: [ val(meta), path(tbi) ]
        filtered_stats         = FILTERMUTECTCALLS.out.stats.collect()                 // channel: [ val(meta), path(stats) ]

        versions               = ch_versions                                           // channel: [ path(versions.yml) ]
}
