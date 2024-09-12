//
// VCF VARIANT ANNOTATION

// Include nf-core subworkflows
include { TABIX_BGZIP } from '../../../modules/nf-core/tabix/bgzip/main.nf'
include { VCF_ANNOTATE_SNPEFF } from '../../../subworkflows/nf-core/vcf_annotate_snpeff/main'
include { VCF_ANNOTATE_ENSEMBLVEP} from '../../../subworkflows/nf-core/vcf_annotate_ensemblvep/main'

// Include nf-core modules

workflow VCF_SHORT_VARIANT_ANNOTATION {
  take:
    ch_vcf  // channel: [ meta, path(vcf) ]
    ch_fasta // channel: [ meta, path(fasta) ]
  main:
    ch_versions = Channel.empty()
    ch_annotated_vcf_tbi = Channel.empty()

    def snpeff_db = params.genomes[params.genome].snpeff_db
    def snpeff_cache = file( params.genomes[params.genome].snpeff_cache, checkIfExists: true )

    def vep_species = params.genomes[params.genome].vep_species
    def vep_cache_version = params.genomes[params.genome].vep_cache_version
    def vep_cache_dir = file( params.genomes[params.genome].vep_cache_dir, checkIfExists: true )
    ch_vep_cache_dir = Channel.value( vep_cache_dir )
      .map{ vep_cache -> [ [ id:'vep_cache_dir' ], vep_cache ] } 

    ch_snpeff_cache = Channel.value( snpeff_cache )
          .map{ cache -> [ [ id:'snpeff_cache' ], cache ] }

    for ( tool in params.vcf_short_variant_annotation.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "snpeff" ) {
        ch_snpeff_vcf = ch_vcf.map{ meta, vcf -> 
          meta = meta + [ann: "snpeff" ]
          meta.id = meta.id+".snpeff"
          [ meta, vcf ]
        }
        VCF_ANNOTATE_SNPEFF( ch_snpeff_vcf, snpeff_db, ch_snpeff_cache )
        ch_versions = ch_versions.mix( VCF_ANNOTATE_SNPEFF.out.versions )
        ch_annotated_vcf_tbi = ch_annotated_vcf_tbi.mix( VCF_ANNOTATE_SNPEFF.out.vcf_tbi )
      
        known_tool = true
      }

      if (tool == "vep" ) {
        // unzip input to avoid broken pipe error: 
        TABIX_BGZIP( ch_vcf )
        ch_vep_vcf = TABIX_BGZIP.out.output
        ch_vep_vcf = ch_vep_vcf.map{ meta, vcf -> 
          meta = meta + [ann: "vep" ]
          meta.id = meta.id+".vep"
          [ meta, vcf, [] ]
        }

        VCF_ANNOTATE_ENSEMBLVEP( ch_vep_vcf, ch_fasta, params.genome, vep_species, vep_cache_version, ch_vep_cache_dir, [] )
        ch_versions = ch_versions.mix( VCF_ANNOTATE_ENSEMBLVEP.out.versions )
        ch_annotated_vcf_tbi = ch_annotated_vcf_tbi.mix( VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi )

        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a vcf variant annotation tool, this tool is not build in (yet).")
      }
    }
  emit:
    versions = ch_versions // channel: [ versions.yml ]
    vcf_tbi = ch_annotated_vcf_tbi // channel: [ meta, vcf, tbi ]
}


