//
// BAM GERMLINE SHORT VARIANT DISCOVERY

// Include local modules
include { GRIDSS } from '../../../modules/local/gridss/gridss/main'
include { GRIDSS_ANNOTATEINSERTEDSEQUENCE } from '../../../modules/local/gridss/annotateinsertedsequence/main'
include { HMFTOOLS_GRIPSS } from '../../../modules/local/hmftools/gripss/main'
include { HMFTOOLS_AMBER } from '../../../modules/local/hmftools/amber/main'
include { HMFTOOLS_COBALT } from '../../../modules/local/hmftools/cobalt/main'
include { HMFTOOLS_PURPLE } from '../../../modules/local/hmftools/purple/main'
include { HMFTOOLS_LINX } from '../../../modules/local/hmftools/linx/linx/main'
include { HMFTOOLS_LINX_SVVISUALISER } from '../../../modules/local/hmftools/linx/svvisualiser/main'

// Include nf-core modules

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX {
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
    def tumor_only_diploid_bed = file( params.genomes[params.genome].cobalt_tumor_only_diploid_bed, checkIfExists: true )

    def driver_gene_panel = file( params.genomes[params.genome].purple_driver_gene_panel, checkIfExists: true )
    def ensembl_data_dir = file( params.genomes[params.genome].purple_ensembl_data_dir, checkIfExists: true )
    def somatic_hotspots = file( params.genomes[params.genome].purple_somatic_hotspots, checkIfExists: true )

    def known_fusion_file = file( params.genomes[params.genome].linx_known_fusion_file, checkIfExists: true )

    ch_bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }

    ch_viralreference_index = Channel.value( file("${params.genomes[params.genome].gridss_viralreference}/*") )
              .map{ viralreference_index -> [ [ id:'index' ], viralreference_index ] }

    GRIDSS( ch_bams, ch_bwa_index, fai, fasta_dict, blacklist, gridss_properties )
    ch_versions = ch_versions.mix( GRIDSS.out.versions )

    ch_gridss_driver_vcfs = GRIDSS.out.gridss_driver_vcf
    GRIDSS_ANNOTATEINSERTEDSEQUENCE( ch_gridss_driver_vcfs, ch_viralreference_index, repeatmaskerbed )
    ch_versions = ch_versions.mix( GRIDSS_ANNOTATEINSERTEDSEQUENCE.out.versions )

    ch_gridss_unfiltered_vcfs = GRIDSS_ANNOTATEINSERTEDSEQUENCE.out.gridss_unfiltered_vcf
    HMFTOOLS_GRIPSS( ch_gridss_unfiltered_vcfs, ch_bwa_index, fai, known_hotspot_file, pon_sgl_file, pon_sv_file, repeat_mask_file )
    ch_versions = ch_versions.mix( HMFTOOLS_GRIPSS.out.versions )

    HMFTOOLS_AMBER( ch_bams, fasta, bafsnps )
    ch_versions = ch_versions.mix( HMFTOOLS_AMBER.out.versions )

    HMFTOOLS_COBALT( ch_bams, fasta, gc_profile, tumor_only_diploid_bed )
    ch_versions = ch_versions.mix( HMFTOOLS_COBALT.out.versions )

    amber = HMFTOOLS_AMBER.out.amber_baf_pcf
      .map{ meta, amber_file ->
        amber_dir = amber_file.getParent()
        [ meta, amber_dir ]
      }
    
    cobalt = HMFTOOLS_COBALT.out.cobalt_gc_median_tsv
      .map{ meta, cobalt_file ->
        cobalt_dir = cobalt_file.getParent()
        [ meta, cobalt_dir ]
      }
    
    ch_gripss_vcf = HMFTOOLS_GRIPSS.out.gripss_vcf
      .join( HMFTOOLS_GRIPSS.out.gripss_vcf_tbi )
    
    ch_gripss_filtered_vcf = HMFTOOLS_GRIPSS.out.gripss_filtered_vcf
      .join( HMFTOOLS_GRIPSS.out.gripss_filtered_vcf_tbi )

    ch_purple = ch_gripss_vcf
      .join( ch_gripss_filtered_vcf )
      .join( amber )
      .join( cobalt )

    HMFTOOLS_PURPLE( ch_purple, fasta, fai, fasta_dict, gc_profile, driver_gene_panel, somatic_hotspots, ensembl_data_dir )
    ch_versions = ch_versions.mix( HMFTOOLS_PURPLE.out.versions )
    
    ch_linx = HMFTOOLS_PURPLE.out.purple_sv_vcf
      .join( HMFTOOLS_PURPLE.out.purple_sv_vcf_tbi )
      .map{ meta, purple_vcf, purple_vcf_tbi ->
        purple_dir = purple_vcf.getParent()
        [ meta, purple_vcf, purple_vcf_tbi, purple_dir]
      }

    HMFTOOLS_LINX( ch_linx, ensembl_data_dir, driver_gene_panel, known_fusion_file )
    ch_versions = ch_versions.mix( HMFTOOLS_LINX.out.versions )

    linx = HMFTOOLS_LINX.out.linx_vis_copy_number_tsv
      .map{ meta, linx_file -> 
        linx_dir = linx_file.getParent()
        [ meta, linx_dir ]
      }

    HMFTOOLS_LINX_SVVISUALISER( linx, ensembl_data_dir )
    ch_versions = ch_versions.mix( HMFTOOLS_LINX_SVVISUALISER.out.versions )

  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


