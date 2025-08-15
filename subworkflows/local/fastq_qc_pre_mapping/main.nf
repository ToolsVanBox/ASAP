//
// PRE MAPPING QC

// Include nf-core modules
include { FASTQC } from '../../../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../../../modules/nf-core/multiqc/main.nf'

workflow FASTQ_QC_PRE_MAPPING {
  take:
    ch_fastq // channel: [ meta, path(fastq) ]
  main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()
    
    def run_fastqc = false 
    
    for ( tool in params.fastq_qc_pre_mapping.tool ) {
      tool = tool.toLowerCase()      
      def known_tool = false
      
      // Run FastQC 
      if ( tool  == "fastqc" ) {

          FASTQC( ch_fastq )
          ch_reports = ch_reports.mix( FASTQC.out.zip.collect{ meta, logs -> logs } )
          ch_versions = ch_versions.mix( FASTQC.out.versions.first() )

          run_fastqc = true          
          known_tool = true
      }
      
      // Run MultiQC
      if ( tool == "multiqc" ) {    
        if ( ! run_fastqc ) {
          error("MultiQC cannot be runned in the pre mapping QC step without FastQC, please add FastQC to params.fastq_qc_pre_mapping.tool in ${projectDir}/configs/subworkflows/local/fastq_qc_pre_mapping")
        } 

        CUSTOM_DUMPSOFTWAREVERSIONS( ch_versions.unique().collectFile(name: 'collated_versions.yml') )
        
        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix( CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect() )
        multiqc_files = multiqc_files.mix( FASTQC.out.zip.collect{it[1]}.ifEmpty([]) )

        multiqc_custom_config = Channel.empty()
        multiqc_logo = Channel.empty()

        multiqc_config = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )

        MULTIQC( multiqc_files.collect(), multiqc_config.toList(), multiqc_custom_config.toList(), multiqc_logo.toList() )

        multiqc_report = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix( MULTIQC.out.versions)
        
        known_tool = true    
    }

    if ( ! known_tool ) {
      println ("WARNING: Skip ${tool}, because it's not known as a QC pre mapping tool, this tool is not build in (yet).")
    }
  }

  emit:
    fastqc = ch_reports // channel: [ fastqc_reports ]
    multiqc = multiqc_report // channel: [ multiqc_report ]
    versions = ch_versions // channel: [ versions.yml ]
        
}
