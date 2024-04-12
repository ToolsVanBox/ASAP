//
// BAM GERMLINE SHORT VARIANT DISCOVERY WITH GATK4

// Include nf-core modules
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'                                                                      
include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'                                                                                                     

// Include local subworkflows
include { BAM_MERGE } from '../../../subworkflows/local/bam_merge/main.nf'

workflow BAM_BASE_RECALIBRATION_GATK4 {
  take:
    ch_bam_interval   // channel: [ val(meta), path(bam), path(bai), path(interval) ]
    ch_fasta   // channel: [ val(meta), path(fasta) ]
    ch_fai   // channel: [ val(meta), path(fai) ]
    ch_dict   // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_bqsr_bam = Channel.empty()
    ch_bqsr_bai = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    fai = ch_fai.map{ meta, fai -> [ fai ] }
    dict = ch_dict.map{ meta, dict -> [ dict ] }

    known_sites = []
    known_sites_tbi = []
    
    for ( known_site_file in params.genomes[params.genome].known_sites ) {
      def known_site_tbi_file = known_site_file+".tbi"
      def file_exists = file(known_site_file).exists()
      def tbi_exists = file(known_site_tbi_file).exists()
      if (! (file_exists | tbi_exists) ) {
        println ">>> WARNING: Known sites file ${known_site_tbi_file} and/or ${known_site_tbi_file} does not exist and will not be included"
      }
      else {
        known_sites.add( known_site_file )
        known_sites_tbi.add( known_site_tbi_file )
      }
    }
      
    GATK4_BASERECALIBRATOR( ch_bam_interval, fasta, fai, dict, known_sites, known_sites_tbi )
    ch_versions = ch_versions.mix( GATK4_BASERECALIBRATOR.out.versions )  

    bqsr_table = GATK4_BASERECALIBRATOR.out.table
    ch_bam_interval_bqsrtable = ch_bam_interval.join( bqsr_table )
      .map{ meta, bam, bai, interval, bqsrtable -> 
        meta = meta + [ bqsr: "gatk" ]
        meta.id = meta.id+".bqsr"
        [ meta, bam, bai, bqsrtable, interval ]
      }

    GATK4_APPLYBQSR( ch_bam_interval_bqsrtable, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_APPLYBQSR.out.versions )  

    ch_bqsr_interval_bams = GATK4_APPLYBQSR.out.bam
      .map{ meta, bam -> 
        meta = meta - meta.subMap("interval")
        meta.id = meta.id.replaceAll(/.\d+.bqsr$/,".bqsr")
        [ meta, bam ]
      }
      .groupTuple()


    BAM_MERGE( ch_bqsr_interval_bams, ch_fasta, ch_fai )
    ch_versions = ch_versions.mix( BAM_MERGE.out.versions )  

    ch_bqsr_bam = ch_bqsr_bam.mix( BAM_MERGE.out.bam )
    ch_bqsr_bai = ch_bqsr_bai.mix( BAM_MERGE.out.bai )

  emit:
    bam = ch_bqsr_bam  // channel: [ meta, bam ]
    bai = ch_bqsr_bai  // channel: [ meta, bai ]
    versions = ch_versions // channel: [ versions.yml ]

}


