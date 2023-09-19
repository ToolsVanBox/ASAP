//
// BAM GERMLINE SHORT VARIANT DISCOVERY

// Include local subworkflows
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX } from '../../../subworkflows/local/bam_germline_structural_variant_discovery_gridss_purple_linx/main'

// Include nf-core modules

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY {
  take:
    ch_bams  // channel: [ meta, path(bam), path(bai) ]
  main:
    ch_versions = Channel.empty()

    BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_GRIDSS_PURPLE_LINX( ch_bams )
    // def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    // def fai = file( fasta.toString()+".fai", checkIfExists: true )
    
    // dbsnp = []
    // dbsnp_tbi = []

    // for ( dbsnp_file in params.genomes[params.genome].dbsnp ) {
    //   def dbsnp_tbi_file = dbsnp_file+".tbi"
    //   def file_exists = file(dbsnp_file).exists()
    //   def tbi_exists = file(dbsnp_tbi_file).exists()
    //   if (! (file_exists | tbi_exists) ) {
    //     println ">>> WARNING: DBSNP file ${dbsnp_file} and/or ${dbsnp_tbi_file} does not exist and will not be included"
    //   }
    //   else {
    //     dbsnp.add( dbsnp_file )
    //     dbsnp_tbi.add( dbsnp_tbi_file )
    //   }
    // }

    // for ( tool in params.bam_germline_copy_number_discovery.tool ) {
    //   tool = tool.toLowerCase()      
    //   known_tool = false    

    //   if ( tool == "controlfreec" ) {
        
    //     BAM_GERMLINE_COPY_NUMBER_DISCOVERY_CONTROLFREEC( ch_bams, fasta, fai, dbsnp, dbsnp_tbi )
    //     ch_versions = ch_versions.mix( BAM_GERMLINE_COPY_NUMBER_DISCOVERY_CONTROLFREEC.out.versions )
      
    //     known_tool = true
    //   }

    //   if ( ! known_tool ) {
    //    println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
    //   }
    // }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


