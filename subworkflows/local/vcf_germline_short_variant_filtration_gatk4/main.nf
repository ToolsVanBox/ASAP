//
// VCF VARIANT FILTRATION WITH GATK4

// Include nf-core modules
include { GATK4_VARIANTFILTRATION } from '../../../modules/nf-core/gatk4/variantfiltration/main'

workflow VCF_GERMLINE_SHORT_VARIANT_FILTRATION_GATK4 {
  take:
    ch_vcf_tbi  // channel: [ meta, path(vcf), path(tbi) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    
  main:
    ch_versions = Channel.empty()
    ch_vcfs = Channel.empty()
    ch_tbi = Channel.empty()

    GATK4_VARIANTFILTRATION( ch_vcf_tbi, ch_fasta, ch_fai, ch_dict )
    ch_versions = ch_versions.mix( GATK4_VARIANTFILTRATION.out.versions )

    ch_vcfs = ch_vcfs.mix( GATK4_VARIANTFILTRATION.out.vcf )
    ch_tbi = ch_tbi.mix( GATK4_VARIANTFILTRATION.out.tbi )

  emit:
    vcf = ch_vcfs // channel: [ meta, vcf ]
    tbi = ch_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


