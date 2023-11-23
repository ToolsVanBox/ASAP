//
// VCF VARIANT ANNOTATION

// Include nf-core subworkflows
include { VCF_ANNOTATE_SNPEFF } from '../../../subworkflows/nf-core/vcf_annotate_snpeff/main'

// Include nf-core modules

workflow VCF_VARIANT_ANNOTATION {
  take:
    ch_vcf  // channel: [ meta, path(vcf) ]
  main:
    ch_versions = Channel.empty()

    def snpeff_db = params.genomes[params.genome].snpeff_db
    def snpeff_cache = file( params.genomes[params.genome].snpeff_cache, checkIfExists: true )

    ch_snpeff_cache = Channel.value( snpeff_cache )
          .map{ cache -> [ [ id:'snpeff_cache' ], cache ] }

    for ( tool in params.vcf_variant_annotation.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "snpeff" ) {
        
        VCF_ANNOTATE_SNPEFF( ch_vcf, snpeff_db, ch_snpeff_cache )
        ch_versions = ch_versions.mix( VCF_ANNOTATE_SNPEFF.out.versions )
      
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf variant annotation tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
}


