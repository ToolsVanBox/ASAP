//
// BAM TUMOR ONLY STRUCTURAL VARIANT DISCOVERY WITH GRIDSS PURPLE LINX

// Include local modules
include { GRIDSS as GRIDSS_TUMORONLY } from '../../../modules/local/gridss/gridss/main'
include { GRIDSS_ANNOTATEINSERTEDSEQUENCE as GRIDSS_ANNOTATEINSERTEDSEQUENCE_TUMORONLY } from '../../../modules/local/gridss/annotateinsertedsequence/main'
include { HMFTOOLS_GRIPSS_TUMORONLY } from '../../../modules/local/hmftools/gripss/tumoronly/main'
include { HMFTOOLS_AMBER_TUMORONLY } from '../../../modules/local/hmftools/amber/tumoronly/main'
include { HMFTOOLS_COBALT_TUMORONLY } from '../../../modules/local/hmftools/cobalt/tumoronly/main'
include { HMFTOOLS_PURPLE_TUMORONLY } from '../../../modules/local/hmftools/purple/tumoronly/main'
include { HMFTOOLS_LINX_TUMORONLY } from '../../../modules/local/hmftools/linx/tumoronly/main'
include { HMFTOOLS_LINX_SVVISUALISER as HMFTOOLS_LINX_SVVISUALISER_TUMORONLY } from '../../../modules/local/hmftools/linx/svvisualiser/main'

// Include nf-core modules

workflow BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX {
  take:
    ch_bams  // channel: [ meta, path(bam), path(bai) ]
  main:
    ch_versions = Channel.empty()

    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( fasta.toString()+".fai", checkIfExists: true )
    def bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
    def fasta_dict = file( fasta.toString()+".dict", checkIfExists: true )

    def blacklist = file( params.genomes[params.genome].gridss_blacklist, checkIfExists: true )
    def gridss_properties = file( params.genomes[params.genome].gridss_properties, checkIfExists: true )
    def repeatmaskerbed = file( params.genomes[params.genome].gridss_repeatmaskerbed, checkIfExists: true )
    
    def known_hotspot_file = file( params.genomes[params.genome].gripss_known_hotspot_file, checkIfExists: true )
    def pon_sgl_file = file( params.genomes[params.genome].gripss_pon_sgl_file, checkIfExists: true )
    def pon_sv_file = file( params.genomes[params.genome].gripss_pon_sv_file, checkIfExists: true )
    def repeat_mask_file = file( params.genomes[params.genome].gripss_repeat_mask_file, checkIfExists: true )

    def bafsnps = file( params.genomes[params.genome].amber_bafsnps, checkIfExists: true )

    def gc_profile = file( params.genomes[params.genome].cobalt_gc_profile, checkIfExists: true )
    def tumoronly_diploid_bed = file( params.genomes[params.genome].cobalt_tumoronly_diploid_bed, checkIfExists: true )

    def driver_gene_panel = file( params.genomes[params.genome].purple_driver_gene_panel, checkIfExists: true )
    def ensembl_data_dir = file( params.genomes[params.genome].purple_ensembl_data_dir, checkIfExists: true )
    def somatic_hotspots = file( params.genomes[params.genome].purple_somatic_hotspots, checkIfExists: true )

    def known_fusion_file = file( params.genomes[params.genome].linx_known_fusion_file, checkIfExists: true )

    ch_bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }

    ch_viralreference_index = Channel.value( file("${params.genomes[params.genome].gridss_viralreference}/*") )
              .map{ viralreference_index -> [ [ id:'index' ], viralreference_index ] }

    ch_gridss_bam = ch_bams.map{ meta, bam, bai -> [ [ id: meta.run_id ], bam ] }.groupTuple()
    ch_gridss_bai = ch_bams.map{ meta, bam, bai -> [ [ id: meta.run_id ], bai ] }.groupTuple()
    ch_gridss_labels = ch_bams.map{ meta, bam, bai -> [ [ id: meta.run_id ], meta.sample_id ] }.groupTuple()
    
    GRIDSS_TUMORONLY( ch_gridss_bam, ch_gridss_bai, ch_gridss_labels, ch_bwa_index, fai, fasta_dict, blacklist, gridss_properties )
    ch_versions = ch_versions.mix( GRIDSS_TUMORONLY.out.versions )

    ch_gridss_driver_vcfs = GRIDSS_TUMORONLY.out.gridss_driver_vcf
    GRIDSS_ANNOTATEINSERTEDSEQUENCE_TUMORONLY( ch_gridss_driver_vcfs, ch_viralreference_index, repeatmaskerbed )
    ch_versions = ch_versions.mix( GRIDSS_ANNOTATEINSERTEDSEQUENCE_TUMORONLY.out.versions )

    ch_gridss_unfiltered_vcfs = GRIDSS_ANNOTATEINSERTEDSEQUENCE_TUMORONLY.out.gridss_unfiltered_vcf

    ch_bams_sample_type = ch_bams.branch{
        normal: it[0].sample_type == "normal"
        tumor: it[0].sample_type == "tumor"
    }

    ch_bams_tumor_normal = ch_bams_sample_type.tumor
            .combine( ch_bams_sample_type.normal.ifEmpty( [ [], null, null ] ) )
    
    ch_bams_tumor = ch_bams_tumor_normal
            .map{ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
                if ( normal_bam == null ) { error("Need a tumor sample for tumor only calling") }
                [ [ id: tumor_meta.id, run_id: tumor_meta.run_id ], tumor_bam, tumor_bai ]
            }       

    ch_gripss = ch_gridss_unfiltered_vcfs.combine( ch_bams_tumor )
      .map{ meta, gridss_unfiltered_vcf, meta2, tumor_bam, tumor_bai ->
        [ meta2, gridss_unfiltered_vcf ]
       }
    
    HMFTOOLS_GRIPSS_TUMORONLY( ch_gripss, ch_bwa_index, fai, known_hotspot_file, pon_sgl_file, pon_sv_file, repeat_mask_file )
    ch_versions = ch_versions.mix( HMFTOOLS_GRIPSS_TUMORONLY.out.versions )

    HMFTOOLS_AMBER_TUMORONLY( ch_bams_tumor, fasta, bafsnps )
    ch_versions = ch_versions.mix( HMFTOOLS_AMBER_TUMORONLY.out.versions )

    HMFTOOLS_COBALT_TUMORONLY( ch_bams_tumor, fasta, gc_profile, tumoronly_diploid_bed )
    ch_versions = ch_versions.mix( HMFTOOLS_COBALT_TUMORONLY.out.versions )

    amber_dir = HMFTOOLS_AMBER_TUMORONLY.out.amber_baf_pcf
      .map{ meta, amber_file ->
        amber_dir = amber_file.getParent()
        [ meta, amber_dir ]
      }
    
    cobalt_dir = HMFTOOLS_COBALT_TUMORONLY.out.cobalt_ratio_tsv_gz
      .map{ meta, cobalt_file ->
        cobalt_dir = cobalt_file.getParent()
        [ meta, cobalt_dir ]
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

  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


