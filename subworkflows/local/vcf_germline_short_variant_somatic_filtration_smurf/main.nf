//
// VCF SOMATIC FILTRATION SMURF

// Include local modules

include { SMURF } from '../../../modules/local/SMuRF/main.nf'

// Include nf-core modules
include { TABIX_BGZIP } from '../../../modules/nf-core/tabix/bgzip/main.nf'
include { SNPSIFT_SPLIT } from '../../../modules/nf-core/snpsift/split/main.nf'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_SMURF_FILTERED_VCFS } from '../../../modules/nf-core/snpsift/split/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_SMURF_VCFS } from '../../../modules/nf-core/snpsift/split/main.nf'

workflow VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION_SMURF {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_bam_bai // channel: [ meta, path(bam), path(bai) ]
    
  main:
    def config = file( params.genomes[params.genome].smurf_config, checkIfExists: true )
     
    ch_input = ch_vcf_tbi
      .map{ meta, vcf, tbi ->
        meta = meta + [ split: true ]
        [ meta, vcf ]
      }

    TABIX_BGZIP( ch_input )

    ch_vcf = TABIX_BGZIP.out.output

    SNPSIFT_SPLIT( ch_vcf )
    
    ch_bgziptabix = SNPSIFT_SPLIT.out.out_vcfs
        .map{ meta, vcf_file ->
            [ vcf_file ]
        }
        .flatten()
        .map{ vcf_file -> 
            [ [ id: vcf_file.getBaseName(), split: true ], vcf_file ]
        }
          
    TABIX_BGZIPTABIX( ch_bgziptabix )

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

    // ch_smurf = TABIX_BGZIPTABIX.out.gz_tbi
    //   .combine( 
    //     ch_bam_bai
    //       .map{ meta, bam, bai ->
    //         if (meta.sample_type == "normal") {
    //           // meta2 = meta - meta.subMap('sample','sample_type','id')
    //           meta2 = [:]
    //           meta2.run_id = meta.run_id
    //           [ meta2, meta.sample ]
    //         }
    //       }
    //       .groupTuple()
    //       .combine( 
    //         ch_bam_bai
    //           .map{ meta, bam, bai ->
    //             // meta = meta - meta.subMap('sample','sample_type','id')
    //             meta2 = [:]
    //             meta2.run_id = meta.run_id
    //             [meta2, bam, bai]
    //           }
    //           .groupTuple()
    //       )
    //       .map{ meta, bulk_names2, meta2, bam, bai ->
    //         meta = meta + [ bulk_names: bulk_names2 ]
    //         [ meta, bam, bai ]
    //       }
    //   )
    //   .map{ meta, vcf_file, tbi, meta2, bam, bai ->
    //       meta = meta + [ bulk_names: meta2.bulk_names ]
    //       [ meta, vcf_file, tbi, bam, bai ]
    //   }

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

      SNPSIFT_JOIN_SMURF_FILTERED_VCFS( ch_filtered_vcfs )

      SNPSIFT_JOIN_SMURF_VCFS( ch_smurf_vcfs )

  // emit:
    
}


