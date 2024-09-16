//
// Run GATK mutect2 in tumor normal mode
//

include { GATK4_MUTECT2                   as MUTECT2 }                   from '../../../modules/local/gatk4/mutect2/main'
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

        ch_merge_mutect2_f1r2 = MUTECT2.out.f1r2
            .map{ meta, f1r2_file -> 
                meta = meta - meta.subMap("id","interval","normal_ids")
                meta = meta + [ id: meta.run_id ]
                [ meta, f1r2_file ]
            }
            .groupTuple()

        MERGEMUTECTSTATS(ch_merge_mutect2_stats)   
        ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)

    emit:
        vcf            = MERGE_MUTECT2.out.vcf.collect()                       // channel: [ val(meta), path(vcf) ]
        tbi            = MERGE_MUTECT2.out.tbi.collect()                       // channel: [ val(meta), path(tbi) ]
        stats          = MERGEMUTECTSTATS.out.stats.collect()                  // channel: [ val(meta), path(stats) ]
        f1r2           = ch_merge_mutect2_f1r2                            // channel: [ val(meta), path(f1r2) ]

        versions       = ch_versions                                           // channel: [ path(versions.yml) ]
}
