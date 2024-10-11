//
// VCF SOMATIC FILTRATION SMURF

// Include local modules
include { SMURF } from '../../../modules/local/SMuRF/main.nf'

// Include nf-core modules
include { TABIX_BGZIP } from '../../../modules/nf-core/tabix/bgzip/main.nf'
include { SNPSIFT_SPLIT } from '../../../modules/nf-core/snpsift/split/main.nf'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_TABIX  as TABIX_SMURF } from '../../../modules/nf-core/tabix/tabix/main.nf'
include { TABIX_TABIX  as TABIX_SMURF_FILTERED } from '../../../modules/nf-core/tabix/tabix/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_SMURF_FILTERED_VCFS } from '../../../modules/nf-core/snpsift/split/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_SMURF_VCFS } from '../../../modules/nf-core/snpsift/split/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SMURF_SORT} from '../../../modules/nf-core/bcftools/sort/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SMURF_FILTERED_SORT } from '../../../modules/nf-core/bcftools/sort/main.nf'

workflow VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ]
    
  main:
    ch_versions = Channel.empty()
    
    def config = file( params.genomes[params.genome].smurf_config, checkIfExists: true )
     
    ch_input = ch_vcf_tbi
      .map{ meta, vcf, tbi ->
        meta = meta + [ split: true ]
        [ meta, vcf ]
      }

    TABIX_BGZIP( ch_input )
    ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )
    ch_vcf = TABIX_BGZIP.out.output

    SNPSIFT_SPLIT( ch_vcf )
    ch_versions = ch_versions.mix( SNPSIFT_SPLIT.out.versions )
    
    ch_bgziptabix = SNPSIFT_SPLIT.out.out_vcfs
        .map{ meta, vcf_file ->
            [ vcf_file ]
        }
        .flatten()
        .map{ vcf_file -> 
            [ [ id: vcf_file.getBaseName(), split: true ], vcf_file ]
        }
          
    TABIX_BGZIPTABIX( ch_bgziptabix )
    ch_versions = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    ch_smurf = TABIX_BGZIPTABIX.out.gz_tbi
      .combine( 
        ch_bam_bai
          .map{ meta, bam, bai ->
            if (meta.sample_type == "normal") {
              // meta2 = meta - meta.subMap('sample','sample_type','id')
              meta2 = [:]
              meta2.run_id = meta.run_id
              [ meta2, meta.sample ]
            } else {
              meta2 = [:]
              meta2.run_id = meta.run_id
              [ meta2, [] ]
            }
          }
          .groupTuple()
          .combine( 
            ch_bam_bai
              .map{ meta, bam, bai ->
                // meta = meta - meta.subMap('sample','sample_type','id')
                meta2 = [:]
                meta2.run_id = meta.run_id
                [meta2, bam, bai]
              }
              .groupTuple()
          )
          .map{ meta, bulk_names2, meta2, bam, bai ->
            meta = meta + [ bulk_names: bulk_names2 ]
            [ meta, bam, bai ]
          }
      )
      .map{ meta, vcf_file, tbi, meta2, bam, bai ->
          meta = meta + [ bulk_names: meta2.bulk_names ]
          [ meta, vcf_file, tbi, bam, bai ]
      }

      SMURF( ch_smurf, config )

      ch_filtered_vcfs = SMURF.out.smurf_filtered_vcf
        .map{ meta, vcf_file -> 
            [ [ id: file(meta.id).getBaseName()+".SMuRF.filtered" ], vcf_file ]
        }
        .groupTuple()
    
      ch_smurf_vcfs = SMURF.out.smurf_vcf
          .map{ meta, vcf_file -> 
              [ [ id: file(meta.id).getBaseName()+".SMuRF" ], vcf_file ]
          }
          .groupTuple()
      
      // Join, sort and tabix for the smurf filtered vcf
      SNPSIFT_JOIN_SMURF_FILTERED_VCFS( ch_filtered_vcfs )
      BCFTOOLS_SMURF_FILTERED_SORT( SNPSIFT_JOIN_SMURF_FILTERED_VCFS.out.out_vcfs )
      ch_versions = ch_versions.mix( BCFTOOLS_SMURF_FILTERED_SORT.out.versions )
      TABIX_SMURF_FILTERED( BCFTOOLS_SMURF_FILTERED_SORT.out.vcf )
      ch_versions = ch_versions.mix( TABIX_SMURF_FILTERED.out.versions )

      // Join, sort and tabix for the smurf vcf
      SNPSIFT_JOIN_SMURF_VCFS( ch_smurf_vcfs )
      BCFTOOLS_SMURF_SORT( SNPSIFT_JOIN_SMURF_VCFS.out.out_vcfs )
      ch_versions = ch_versions.mix( BCFTOOLS_SMURF_SORT.out.versions )
      TABIX_SMURF( BCFTOOLS_SMURF_SORT.out.vcf )
      ch_versions = ch_versions.mix( TABIX_SMURF.out.versions )
      
      
      

  emit:
    versions = ch_versions                                // channel: [ versions.yml ]

    filtered_vcf = BCFTOOLS_SMURF_FILTERED_SORT.out.vcf   // channel: [ meta, filtered_vcf ]
    filtered_tbi = TABIX_SMURF_FILTERED.out.tbi           // channel: [ meta, filtered_tbi ]

    vcf = BCFTOOLS_SMURF_SORT.out.vcf                     // channel: [ meta, vcf ]
    tbi = TABIX_SMURF.out.tbi                             // channel: [ meta, tbi ]
    
}


