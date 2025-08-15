//
// BAM TUMOR-ONLY SHORT VARIANT DISCOVERY

// Include local subworkflows
include { BAM_TUMORONLY_COPY_NUMBER_DISCOVERY_CONTROLFREEC } from '../../../subworkflows/local/bam_tumoronly_copy_number_discovery_controlfreec/main'

// Include nf-core modules

workflow BAM_TUMORONLY_COPY_NUMBER_DISCOVERY {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
  main:
    ch_versions = Channel.empty()
    
    dbsnp = []
    dbsnp_tbi = []

    for ( dbsnp_file in params.genomes[params.genome].dbsnp ) {
      def dbsnp_tbi_file = dbsnp_file+".tbi"
      def file_exists = file(dbsnp_file).exists()
      def tbi_exists = file(dbsnp_tbi_file).exists()
      if (! (file_exists | tbi_exists) ) {
        println ">>> WARNING: DBSNP file ${dbsnp_file} and/or ${dbsnp_tbi_file} does not exist and will not be included"
      }
      else {
        dbsnp.add( dbsnp_file )
        dbsnp_tbi.add( dbsnp_tbi_file )
      }
    }

    for ( tool in params.bam_tumoronly_copy_number_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "controlfreec" ) {
        
        BAM_TUMORONLY_COPY_NUMBER_DISCOVERY_CONTROLFREEC( ch_bam_bai, ch_fasta, dbsnp, dbsnp_tbi )
        ch_versions = ch_versions.mix( BAM_TUMORONLY_COPY_NUMBER_DISCOVERY_CONTROLFREEC.out.versions )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam tumor-only short variant discovery tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


