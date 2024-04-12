//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY WITH GRIDSS PURPLE LINX

// Include local modules
include { HMFTOOLS_GRIPSS_GERMLINE } from '../../../modules/local/hmftools/gripss/germline/main'
include { HMFTOOLS_AMBER_GERMLINE } from '../../../modules/local/hmftools/amber/germline/main'
include { HMFTOOLS_COBALT_GERMLINE } from '../../../modules/local/hmftools/cobalt/germline/main'
include { HMFTOOLS_PURPLE_GERMLINE } from '../../../modules/local/hmftools/purple/germline/main'

include { HMFTOOLS_GRIPSS_SOMATIC } from '../../../modules/local/hmftools/gripss/somatic/main'
include { HMFTOOLS_AMBER_SOMATIC } from '../../../modules/local/hmftools/amber/somatic/main'
include { HMFTOOLS_COBALT_SOMATIC } from '../../../modules/local/hmftools/cobalt/somatic/main'
include { HMFTOOLS_PURPLE_SOMATIC } from '../../../modules/local/hmftools/purple/somatic/main'
include { HMFTOOLS_LINX_SOMATIC } from '../../../modules/local/hmftools/linx/somatic/main'
include { HMFTOOLS_LINX_SVVISUALISER as HMFTOOLS_LINX_SVVISUALISER_SOMATIC } from '../../../modules/local/hmftools/linx/svvisualiser/main'

include { HMFTOOLS_GRIPSS_TUMORONLY } from '../../../modules/local/hmftools/gripss/tumoronly/main'
include { HMFTOOLS_AMBER_TUMORONLY } from '../../../modules/local/hmftools/amber/tumoronly/main'
include { HMFTOOLS_COBALT_TUMORONLY } from '../../../modules/local/hmftools/cobalt/tumoronly/main'
include { HMFTOOLS_PURPLE_TUMORONLY } from '../../../modules/local/hmftools/purple/tumoronly/main'
include { HMFTOOLS_LINX_TUMORONLY } from '../../../modules/local/hmftools/linx/tumoronly/main'
include { HMFTOOLS_LINX_SVVISUALISER as HMFTOOLS_LINX_SVVISUALISER_TUMORONLY } from '../../../modules/local/hmftools/linx/svvisualiser/main'

// Include nf-core modules

workflow VCF_STRUCTURAL_VARIANT_FILTRATION_GRIDSS_PURPLE_LINX {
  take:
    ch_vcf  // channel: [ meta, path(vcf) ]
    ch_bam_bai_normal // channel: [ val(meta), path(bam), path(bai) ]
    ch_bam_bai_tumor // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_fasta_dict // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty()

    def fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    def fai = ch_fai.map{ meta, fai -> [ fai ] }
    def fasta_dict = ch_fasta_dict.map{ meta, fasta_dict -> [ fasta_dict ] }

    def known_hotspot_file = file( params.genomes[params.genome].gripss_known_hotspot_file, checkIfExists: true )
    def pon_sgl_file = file( params.genomes[params.genome].gripss_pon_sgl_file, checkIfExists: true )
    def pon_sv_file = file( params.genomes[params.genome].gripss_pon_sv_file, checkIfExists: true )
    def repeat_mask_file = file( params.genomes[params.genome].gripss_repeat_mask_file, checkIfExists: true )
    def bafsnps = file( params.genomes[params.genome].amber_bafsnps, checkIfExists: true )
    def gc_profile = file( params.genomes[params.genome].cobalt_gc_profile, checkIfExists: true )

    def driver_gene_panel = file( params.genomes[params.genome].purple_driver_gene_panel, checkIfExists: true )
    def ensembl_data_dir = file( params.genomes[params.genome].purple_ensembl_data_dir, checkIfExists: true )
    def somatic_hotspots = file( params.genomes[params.genome].purple_somatic_hotspots, checkIfExists: true )
    def germline_hotspots = file( params.genomes[params.genome].purple_germline_hotspots, checkIfExists: true )
    def del_freq_file = file( params.genomes[params.genome].purple_del_freq_file, checkIfExists: true )

    def known_fusion_file = file( params.genomes[params.genome].linx_known_fusion_file, checkIfExists: true )

    
    ch_viralreference_index = Channel.value( file("${params.genomes[params.genome].gridss_viralreference}/*") )
              .map{ viralreference_index -> [ [ id:'index' ], viralreference_index ] }

    ch_bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }

    for ( method in params.vcf_structural_variant_filtration.method ) {
      method = method.toLowerCase()      
      known_method = false    

      if ( method == "germline" ) {
        ch_gripss_germline = ch_bam_bai_normal.ifEmpty( [ [], null, null ] )
          .combine( ch_bam_bai_tumor.ifEmpty( [ [], null, null ] ) )
          .map{ normal_meta, normal_bam, normal_bai, tumor_meta, tumor_bam, tumor_bai ->
              if ( normal_bam == null ) { error("Need a normal sample for structural variant germline filtering") }
              if ( tumor_bam == null ) { error("Need a tumor sample for structutal variant germline calling") }
              [ [ id: "${normal_meta.id}_${tumor_meta.id}", run_id: normal_meta.run_id, tumor_sample_id: tumor_meta.id, normal_sample_id: normal_meta.id ] ]
          }
          .combine( ch_vcf )
          .map{ meta, meta2, vcf -> [ meta, vcf ] }
        
        HMFTOOLS_GRIPSS_GERMLINE( ch_gripss_germline, ch_bwa_index, fai, known_hotspot_file, pon_sgl_file, pon_sv_file, repeat_mask_file )
        ch_versions = ch_versions.mix( HMFTOOLS_GRIPSS_GERMLINE.out.versions )

        HMFTOOLS_AMBER_GERMLINE( ch_bam_bai_normal, fasta, bafsnps )
        ch_versions = ch_versions.mix( HMFTOOLS_AMBER_GERMLINE.out.versions )

        HMFTOOLS_COBALT_GERMLINE( ch_bam_bai_normal, fasta, gc_profile )
        ch_versions = ch_versions.mix( HMFTOOLS_COBALT_GERMLINE.out.versions )

        amber_dir = HMFTOOLS_AMBER_GERMLINE.out.amber_baf_tsv_gz
          .map{ meta, amber_file ->
            amber_dir = amber_file.getParent()
            [ amber_dir ]
          }
        
        cobalt_dir = HMFTOOLS_COBALT_GERMLINE.out.cobalt_gc_median_tsv
          .map{ meta, cobalt_file ->
            cobalt_dir = cobalt_file.getParent()
            [ cobalt_dir ]
          }
        
        ch_gripss_vcf = HMFTOOLS_GRIPSS_GERMLINE.out.gripss_vcf
          .join( HMFTOOLS_GRIPSS_GERMLINE.out.gripss_vcf_tbi )
        
        ch_gripss_filtered_vcf = HMFTOOLS_GRIPSS_GERMLINE.out.gripss_filtered_vcf
          .join( HMFTOOLS_GRIPSS_GERMLINE.out.gripss_filtered_vcf_tbi )
        
        ch_purple = ch_gripss_vcf
          .join( ch_gripss_filtered_vcf )
          .combine( amber_dir )
          .combine( cobalt_dir )
          .map{ meta, gripss_vcf, gripss_tbi, gripss_filtered_vcf, gripss_filtered_vcf_tbi, amber, cobalt ->
            meta.id = meta.normal_sample_id
            [ meta, gripss_vcf, gripss_tbi, gripss_filtered_vcf, gripss_filtered_vcf_tbi, amber, cobalt ] 
          }
          
        HMFTOOLS_PURPLE_GERMLINE( ch_purple, fasta, fai, fasta_dict, gc_profile, driver_gene_panel, germline_hotspots, del_freq_file, ensembl_data_dir )
        ch_versions = ch_versions.mix( HMFTOOLS_PURPLE_GERMLINE.out.versions )

        known_method = true
      }

      if ( method == "somatic" ) {
        ch_gripss_somatic = ch_bam_bai_normal.ifEmpty( [ [], null, null ] )
          .combine( ch_bam_bai_tumor.ifEmpty( [ [], null, null ] ) )
          .map{ normal_meta, normal_bam, normal_bai, tumor_meta, tumor_bam, tumor_bai ->
              if ( normal_bam == null ) { error("Need a normal sample for structural variant somatic filtering") }
              if ( tumor_bam == null ) { error("Need a tumor sample for structutal variant somatic calling") }
              [ [ id: "${normal_meta.id}_${tumor_meta.id}", run_id: normal_meta.run_id, tumor_sample_id: tumor_meta.id, normal_sample_id: normal_meta.id ] ]
          }
          .combine( ch_vcf )
          .map{ meta, meta2, vcf -> [ meta, vcf ] }
        
        HMFTOOLS_GRIPSS_SOMATIC( ch_gripss_somatic, ch_bwa_index, fai, known_hotspot_file, pon_sgl_file, pon_sv_file, repeat_mask_file )
        ch_versions = ch_versions.mix( HMFTOOLS_GRIPSS_SOMATIC.out.versions )

        ch_bam_bai_normal_tumor = ch_bam_bai_normal.ifEmpty( [ [], null, null ] )
          .combine( ch_bam_bai_tumor.ifEmpty( [ [], null, null ] ) )
          .map{ normal_meta, normal_bam, normal_bai, tumor_meta, tumor_bam, tumor_bai ->
              if ( normal_bam == null ) { error("Need a normal sample for structural variant somatic filtering") }
              if ( tumor_bam == null ) { error("Need a tumor sample for structutal variant somatic calling") }
              [ [ id: "${normal_meta.id}_${tumor_meta.id}", run_id: normal_meta.run_id, tumor_sample_id: tumor_meta.id, normal_sample_id: normal_meta.id ], normal_bam, normal_bai, tumor_bam, tumor_bai ]
          }

        HMFTOOLS_AMBER_SOMATIC( ch_bam_bai_normal_tumor, fasta, bafsnps )
        ch_versions = ch_versions.mix( HMFTOOLS_AMBER_SOMATIC.out.versions )

        HMFTOOLS_COBALT_SOMATIC( ch_bam_bai_normal_tumor, fasta, gc_profile )
        ch_versions = ch_versions.mix( HMFTOOLS_COBALT_SOMATIC.out.versions )

        amber_dir = HMFTOOLS_AMBER_SOMATIC.out.amber_baf_tsv_gz
          .map{ meta, amber_file ->
            amber_dir = amber_file.getParent()
            [ meta, amber_dir ]
          }
                
        cobalt_dir = HMFTOOLS_COBALT_SOMATIC.out.cobalt_ratio_tsv_gz
          .map{ meta, cobalt_file ->
            cobalt_dir = cobalt_file.getParent()
            [ meta, cobalt_dir ]
          }
        
        ch_gripss_vcf = HMFTOOLS_GRIPSS_SOMATIC.out.gripss_vcf
          .join( HMFTOOLS_GRIPSS_SOMATIC.out.gripss_vcf_tbi )
        
        ch_gripss_filtered_vcf = HMFTOOLS_GRIPSS_SOMATIC.out.gripss_filtered_vcf
          .join( HMFTOOLS_GRIPSS_SOMATIC.out.gripss_filtered_vcf_tbi )
        
        ch_purple = ch_gripss_vcf
          .join( ch_gripss_filtered_vcf )
          .join( amber_dir )
          .join( cobalt_dir )
          
        HMFTOOLS_PURPLE_SOMATIC( ch_purple, fasta, fai, fasta_dict, gc_profile, driver_gene_panel, somatic_hotspots, ensembl_data_dir )
        ch_versions = ch_versions.mix( HMFTOOLS_PURPLE_SOMATIC.out.versions )
        
        ch_linx = HMFTOOLS_PURPLE_SOMATIC.out.purple_sv_vcf
          .join( HMFTOOLS_PURPLE_SOMATIC.out.purple_sv_vcf_tbi )
          .map{ meta, purple_vcf, purple_vcf_tbi ->
            purple_dir = purple_vcf.getParent()
            [ meta, purple_vcf, purple_vcf_tbi, purple_dir ]
          }

        HMFTOOLS_LINX_SOMATIC( ch_linx, ensembl_data_dir, driver_gene_panel, known_fusion_file )
        ch_versions = ch_versions.mix( HMFTOOLS_LINX_SOMATIC.out.versions )

        linx = HMFTOOLS_LINX_SOMATIC.out.linx_vis_copy_number_tsv
          .map{ meta, linx_file -> 
            linx_dir = linx_file.getParent()
            [ meta, linx_dir ]
          }

        HMFTOOLS_LINX_SVVISUALISER_SOMATIC( linx, ensembl_data_dir )
        ch_versions = ch_versions.mix( HMFTOOLS_LINX_SVVISUALISER_SOMATIC.out.versions )

        known_method = true
      }

      if ( method == 'tumoronly' ) {
        def tumoronly_diploid_bed = file( params.genomes[params.genome].cobalt_tumoronly_diploid_bed, checkIfExists: true )
        
        ch_gripss_tumoronly = ch_bam_bai_tumor.ifEmpty( [ [], null, null ] )
          .map{ tumor_meta, tumor_bam, tumor_bai ->
              if ( tumor_bam == null ) { error("Need a tumor sample for structutal variant somatic calling") }
              [ [ id: "${tumor_meta.id}", run_id: tumor_meta.run_id ] ]
          }
          .combine( ch_vcf )
          .map{ meta, meta2, vcf -> [ meta, vcf ] }
        
        HMFTOOLS_GRIPSS_TUMORONLY( ch_gripss_tumoronly, ch_bwa_index, fai, known_hotspot_file, pon_sgl_file, pon_sv_file, repeat_mask_file )
        ch_versions = ch_versions.mix( HMFTOOLS_GRIPSS_TUMORONLY.out.versions )

        HMFTOOLS_AMBER_TUMORONLY( ch_bam_bai_tumor, fasta, bafsnps )
        ch_versions = ch_versions.mix( HMFTOOLS_AMBER_TUMORONLY.out.versions )

        HMFTOOLS_COBALT_TUMORONLY( ch_bam_bai_tumor, fasta, gc_profile, tumoronly_diploid_bed )
        ch_versions = ch_versions.mix( HMFTOOLS_COBALT_TUMORONLY.out.versions )

        amber_dir = HMFTOOLS_AMBER_TUMORONLY.out.amber_baf_pcf
          .map{ meta, amber_file ->
            amber_dir = amber_file.getParent()            
            [ [ id: "${meta.id}", run_id: meta.run_id ], amber_dir ]
          }
        
        cobalt_dir = HMFTOOLS_COBALT_TUMORONLY.out.cobalt_ratio_tsv_gz
          .map{ meta, cobalt_file ->
            cobalt_dir = cobalt_file.getParent()
            [ [ id: "${meta.id}", run_id: meta.run_id ], cobalt_dir ]
          }
        
        ch_gripss_vcf = HMFTOOLS_GRIPSS_TUMORONLY.out.gripss_vcf
          .join( HMFTOOLS_GRIPSS_TUMORONLY.out.gripss_vcf_tbi )
        
        ch_gripss_filtered_vcf = HMFTOOLS_GRIPSS_TUMORONLY.out.gripss_filtered_vcf
          .join( HMFTOOLS_GRIPSS_TUMORONLY.out.gripss_filtered_vcf_tbi )
        
        ch_purple = ch_gripss_vcf
          .join( ch_gripss_filtered_vcf )
          .join( amber_dir )
          .join( cobalt_dir )

        HMFTOOLS_PURPLE_TUMORONLY( ch_purple, fasta, fai, fasta_dict, gc_profile, driver_gene_panel, somatic_hotspots, ensembl_data_dir )
        ch_versions = ch_versions.mix( HMFTOOLS_PURPLE_TUMORONLY.out.versions )
        
        ch_linx = HMFTOOLS_PURPLE_TUMORONLY.out.purple_sv_vcf
          .join( HMFTOOLS_PURPLE_TUMORONLY.out.purple_sv_vcf_tbi )
          .map{ meta, purple_vcf, purple_vcf_tbi ->
            purple_dir = purple_vcf.getParent()
            [ meta, purple_vcf, purple_vcf_tbi, purple_dir]
          }

        HMFTOOLS_LINX_TUMORONLY( ch_linx, ensembl_data_dir, driver_gene_panel, known_fusion_file )
        ch_versions = ch_versions.mix( HMFTOOLS_LINX_TUMORONLY.out.versions )

        linx = HMFTOOLS_LINX_TUMORONLY.out.linx_vis_copy_number_tsv
          .map{ meta, linx_file -> 
            linx_dir = linx_file.getParent()
            [ meta, linx_dir ]
          }

        HMFTOOLS_LINX_SVVISUALISER_TUMORONLY( linx, ensembl_data_dir )
        ch_versions = ch_versions.mix( HMFTOOLS_LINX_SVVISUALISER_TUMORONLY.out.versions )

        known_method = true
      }

      if ( ! known_method ) {
       println ("WARNING: Skip ${method}, because it's not known as a vcf structural variant filtration method, this method is not build in (yet).")
      }
    }



    // def bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
    
    // def blacklist = file( params.genomes[params.genome].gridss_blacklist, checkIfExists: true )
    // def gridss_properties = file( params.genomes[params.genome].gridss_properties, checkIfExists: true )
    // def repeatmaskerbed = file( params.genomes[params.genome].gridss_repeatmaskerbed, checkIfExists: true )
    
    
    // def bafsnps = file( params.genomes[params.genome].amber_bafsnps, checkIfExists: true )

    // def gc_profile = file( params.genomes[params.genome].cobalt_gc_profile, checkIfExists: true )
    // def tumor_only_diploid_bed = file( params.genomes[params.genome].cobalt_tumor_only_diploid_bed, checkIfExists: true )

    // def driver_gene_panel = file( params.genomes[params.genome].purple_driver_gene_panel, checkIfExists: true )
    // def ensembl_data_dir = file( params.genomes[params.genome].purple_ensembl_data_dir, checkIfExists: true )
    // def somatic_hotspots = file( params.genomes[params.genome].purple_somatic_hotspots, checkIfExists: true )
    // def germline_hotspots = file( params.genomes[params.genome].purple_germline_hotspots, checkIfExists: true )
    // def del_freq_file = file( params.genomes[params.genome].purple_del_freq_file, checkIfExists: true )

    // def known_fusion_file = file( params.genomes[params.genome].linx_known_fusion_file, checkIfExists: true )

    
    // ch_viralreference_index = Channel.value( file("${params.genomes[params.genome].gridss_viralreference}/*") )
    //           .map{ viralreference_index -> [ [ id:'index' ], viralreference_index ] }

    // ch_gridss_bam = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], bam ] }.groupTuple()
    // ch_gridss_bai = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], bai ] }.groupTuple()
    // ch_gridss_labels = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], meta.sample_id ] }.groupTuple()
    
    // GRIDSS_GERMLINE( ch_gridss_bam, ch_gridss_bai, ch_gridss_labels, ch_bwa_index, fai, fasta_dict, blacklist, gridss_properties )
    // ch_versions = ch_versions.mix( GRIDSS_GERMLINE.out.versions )

    // ch_gridss_driver_vcfs = GRIDSS_GERMLINE.out.gridss_driver_vcf
    // GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE( ch_gridss_driver_vcfs, ch_viralreference_index, repeatmaskerbed )
    // ch_versions = ch_versions.mix( GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE.out.versions )

    // ch_gridss_unfiltered_vcfs = GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE.out.gridss_unfiltered_vcf
    
    // ch_bam_bai_sample_type = ch_bam_bai.branch{
    //     normal: it[0].sample_type == "normal"
    //     tumor: it[0].sample_type == "tumor"
    // }

    // ch_bam_bai_tumor_normal = ch_bam_bai_sample_type.tumor
    //         .combine( ch_bam_bai_sample_type.normal.ifEmpty( [ [], null, null ] ) )
    
    // ch_bam_bai_normal = ch_bam_bai_tumor_normal
    //         .map{ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
    //             if ( normal_bam == null ) { error("Need a normal sample for germline calling") }
    //             [ [ id: normal_meta.id, run_id: normal_meta.run_id ], normal_bam, normal_bai ]
    //         }
    
    // ch_gripss_germline = ch_gridss_unfiltered_vcfs.combine( ch_bam_bai_normal )
    //   .map{ meta, gridss_unfiltered_vcf, meta2, normal_bam, normal_bai ->
    //     [ meta2, gridss_unfiltered_vcf ]
    //    }

        
    // HMFTOOLS_AMBER_GERMLINE( ch_bam_bai_normal, fasta, bafsnps )
    // ch_versions = ch_versions.mix( HMFTOOLS_AMBER_GERMLINE.out.versions )

    // HMFTOOLS_COBALT_GERMLINE( ch_bam_bai_normal, fasta, gc_profile )
    // ch_versions = ch_versions.mix( HMFTOOLS_COBALT_GERMLINE.out.versions )

    // amber_dir = HMFTOOLS_AMBER_GERMLINE.out.amber_baf_tsv_gz
    //   .map{ meta, amber_file ->
    //     amber_dir = amber_file.getParent()
    //     [ meta, amber_dir ]
    //   }
    
    // cobalt_dir = HMFTOOLS_COBALT_GERMLINE.out.cobalt_gc_median_tsv
    //   .map{ meta, cobalt_file ->
    //     cobalt_dir = cobalt_file.getParent()
    //     [ meta, cobalt_dir ]
    //   }
    
    // ch_gripss_vcf = HMFTOOLS_GRIPSS_GERMLINE.out.gripss_vcf
    //   .join( HMFTOOLS_GRIPSS_GERMLINE.out.gripss_vcf_tbi )
    
    // ch_gripss_filtered_vcf = HMFTOOLS_GRIPSS_GERMLINE.out.gripss_filtered_vcf
    //   .join( HMFTOOLS_GRIPSS_GERMLINE.out.gripss_filtered_vcf_tbi )
    
    // ch_purple = ch_gripss_vcf
    //   .join( ch_gripss_filtered_vcf )
    //   .join( amber_dir )
    //   .join( cobalt_dir )

    // HMFTOOLS_PURPLE_GERMLINE( ch_purple, fasta, fai, fasta_dict, gc_profile, driver_gene_panel, germline_hotspots, del_freq_file, ensembl_data_dir )
    // ch_versions = ch_versions.mix( HMFTOOLS_PURPLE_GERMLINE.out.versions )

  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


