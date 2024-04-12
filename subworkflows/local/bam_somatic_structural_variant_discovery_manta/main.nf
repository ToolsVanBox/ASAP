//
// BAM SOMATIC STRUCTURAL VARIANT DISCOVERY WITH MANTA

// Include nf-core modules
include { MANTA_SOMATIC  } from '../../../modules/nf-core/manta/somatic/main'

// Include nf-core modules

workflow BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY_MANTA {
  take:
    ch_bam_bai_normal_tumor  // channel: [ meta, path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

  main:
    ch_versions = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] } 
    fai = ch_fai.map{ meta, fai -> [ fai ] } 

    MANTA_SOMATIC( ch_bam_bai_normal_tumor, fasta, fai )
    ch_versions = ch_versions.mix( MANTA_SOMATIC.out.versions )
    
    emit:
      versions = ch_versions // channel: [ versions.yml ]
}


