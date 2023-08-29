//
// PRE MAPPING QC

include { FASTQC } from '../../../modules/nf-core/fastqc/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { MULTIQC } from '../../../modules/nf-core/multiqc/main.nf'

workflow FASTQ_QC_PRE_MAPPING {
  take:
    fastq
  main:
    reports = Channel.empty()
    versions = Channel.empty()

    // Convert input tools to UpperCase
    fastqc_qc_pre_mapping_tools = params.fastq_qc_pre_mapping.tool.toString().toUpperCase()

    // Run FastQC 
    if ( fastqc_qc_pre_mapping_tools  == "FASTQC" | fastqc_qc_pre_mapping_tools.contains("FASTQC") ) {
        FASTQC( fastq )
        reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
        versions = versions.mix(FASTQC.out.versions.first())
    }
    
    // Run MultiQC
    if ( fastqc_qc_pre_mapping_tools == "MULTIQC" | fastqc_qc_pre_mapping_tools.contains("MULTIQC") ) {

      if (!(fastqc_qc_pre_mapping_tools == "FASTQC" | fastqc_qc_pre_mapping_tools.contains("FASTQC"))) {
        error("MultiQC cannot be runned in the pre mapping QC step without FastQC, please add FastQC to params.fastq_qc_pre_mapping.tool in ${projectDir}/configs/subworkflows/local/fastq_qc_pre_mapping")
      } 

      CUSTOM_DUMPSOFTWAREVERSIONS( versions.unique().collectFile(name: 'collated_versions.yml') )
      
      multiqc_files = Channel.empty()
      multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
      multiqc_files = multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

      multiqc_custom_config = Channel.empty()
      multiqc_logo = Channel.empty()

      multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

      MULTIQC( multiqc_files.collect(), multiqc_config.toList(), multiqc_custom_config.toList(), multiqc_logo.toList() )

      multiqc_report = MULTIQC.out.report.toList()
      versions = versions.mix(MULTIQC.out.versions)
    }

  emit:
    reports
    versions
}
