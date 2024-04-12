//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY WITH GRIDSS PURPLE LINX

// Include local modules
include { GRIDSS as GRIDSS_GERMLINE } from '../../../modules/local/gridss/gridss/main'
include { GRIDSS_ANNOTATEINSERTEDSEQUENCE as GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE } from '../../../modules/local/gridss/annotateinsertedsequence/main'

// Include nf-core modules

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_fasta_dict // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty()

    def fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    def fai = ch_fai.map{ meta, fai -> [ fai ] }
    def fasta_dict = ch_fasta_dict.map{ meta, fasta_dict -> [ fasta_dict ] }
    
    def bwa_index_folder = new File("${params.genomes[params.genome].bwa}")
    
    def blacklist = file( params.genomes[params.genome].gridss_blacklist, checkIfExists: true )
    def gridss_properties = file( params.genomes[params.genome].gridss_properties, checkIfExists: true )
    def repeatmaskerbed = file( params.genomes[params.genome].gridss_repeatmaskerbed, checkIfExists: true )
  
    ch_bwa_index = Channel.value( file("${params.genomes[params.genome].bwa}/*") )
              .map{ bwa_index -> [ [ id:'index' ], bwa_index ] }

    ch_viralreference_index = Channel.value( file("${params.genomes[params.genome].gridss_viralreference}/*") )
              .map{ viralreference_index -> [ [ id:'index' ], viralreference_index ] }

    ch_gridss_bam = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], bam ] }.groupTuple()
    ch_gridss_bai = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], bai ] }.groupTuple()
    ch_gridss_labels = ch_bam_bai.map{ meta, bam, bai -> [ [ id: meta.run_id ], meta.id ] }.groupTuple()
    
    GRIDSS_GERMLINE( ch_gridss_bam, ch_gridss_bai, ch_gridss_labels, ch_bwa_index, fai, fasta_dict, blacklist, gridss_properties )
    ch_versions = ch_versions.mix( GRIDSS_GERMLINE.out.versions )

    ch_gridss_driver_vcfs = GRIDSS_GERMLINE.out.gridss_driver_vcf
    GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE( ch_gridss_driver_vcfs, ch_viralreference_index, repeatmaskerbed )
    ch_versions = ch_versions.mix( GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE.out.versions )

    ch_gridss_unfiltered_vcfs = GRIDSS_ANNOTATEINSERTEDSEQUENCE_GERMLINE.out.gridss_unfiltered_vcf

  emit:
    vcf = ch_gridss_unfiltered_vcfs // channel: [ meta, unfiltered_vcf ] 
    versions = ch_versions // channel: [ versions.yml ]
}


