//
// BAM GERMLINE STRUCTURAL VARIANT DISCOVERY WITH GRIDSS PURPLE LINX

// Include nf-core modules

include { MANTA_GERMLINE  } from '../../../modules/nf-core/manta/germline/main'

workflow BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY_MANTA {
    take:
      ch_bam_bai       // channel: [ meta, bam, bai ]
    
    main:
      ch_versions = Channel.empty()

      ch_manta = ch_bam_bai
        .map{ meta, bam, bai ->
          [ [id: meta.run_id ], bam, bai ]
        }
        .groupTuple()
        .map { meta, bam, bai -> 
          [ meta, bam, bai, [], [] ]
        }

      def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
      def fai = file( fasta.toString()+".fai", checkIfExists: true )
      
      ch_fasta = Channel.value( fasta )
        .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] } 
      ch_fai = Channel.value( fai )
        .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] } 

      MANTA_GERMLINE( ch_manta, ch_fasta, ch_fai )
      ch_versions = ch_versions.mix( MANTA_GERMLINE.out.versions )

    emit:
      versions = ch_versions // channel: [ versions.yml ]

}