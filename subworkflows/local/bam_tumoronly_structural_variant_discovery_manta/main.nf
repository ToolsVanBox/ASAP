//
// BAM TUMOR ONLY STRUCTURAL VARIANT DISCOVERY WITH MANTA

// Include nf-core modules
include { MANTA_TUMORONLY  } from '../../../modules/nf-core/manta/tumoronly/main'

workflow BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY_MANTA {
  take:
    ch_bams  // channel: [ meta, path(bam), path(bai) ]
  main:
    ch_versions = Channel.empty()

    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( fasta.toString()+".fai", checkIfExists: true )
    
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] } 
    ch_fai = Channel.value( fai )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] } 

    ch_bams_sample_type = ch_bams.branch{
        normal: it[0].sample_type == "normal"
        tumor: it[0].sample_type == "tumor"
    }

    ch_bams_tumor_normal = ch_bams_sample_type.tumor
            .combine( ch_bams_sample_type.normal.ifEmpty( [ [], null, null ] ) )
    
    ch_bams_tumor = ch_bams_tumor_normal
            .map{ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
                if ( tumor_bam == null ) { error("Need a tumor sample for tumoronly calling") }
                [ [ id: tumor_meta.id, run_id: tumor_meta.run_id ], tumor_bam, tumor_bai, [], [] ]
            }       

    MANTA_TUMORONLY( ch_bams_tumor, fasta, fai )
    ch_versions = ch_versions.mix( MANTA_TUMORONLY.out.versions )
    
    emit:
      versions = ch_versions // channel: [ versions.yml ]
}




