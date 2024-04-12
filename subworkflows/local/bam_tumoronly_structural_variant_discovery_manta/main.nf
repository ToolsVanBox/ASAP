//
// BAM TUMOR ONLY STRUCTURAL VARIANT DISCOVERY WITH MANTA

// Include nf-core modules
include { MANTA_TUMORONLY  } from '../../../modules/nf-core/manta/tumoronly/main'

workflow BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_MANTA {
  take:
    ch_bam_bai_tumor  // channel: [ meta, path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

  main:
    ch_versions = Channel.empty()

    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] } 
    fai = ch_fai.map{ meta, fai -> [ fai ] } 

    ch_manta_input = ch_bam_bai_tumor
      .map{ meta, bam, bai -> 
       [ meta, bam, bai, [], [] ]
      }

    MANTA_TUMORONLY( ch_manta_input, fasta, fai )
    ch_versions = ch_versions.mix( MANTA_TUMORONLY.out.versions )
    
    emit:
      versions = ch_versions // channel: [ versions.yml ]
}




