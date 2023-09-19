//
// BAM GERMLINE SHORT VARIANT DISCOVERY

// Include local subworkflows
include { BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER } from '../../../subworkflows/local/bam_germline_short_variant_discovery_gatk4haplotypecaller/main'       
include { BAM_BASE_RECALIBRATION_GATK4 } from '../../../subworkflows/local/bam_base_recalibration_gatk4/main'

// Include nf-core modules
include { GATK4_SPLITINTERVALS } from '../../../modules/nf-core/gatk4/splitintervals/main'                                                                          

workflow BAM_GERMLINE_SHORT_VARIANT_DISCOVERY {
  take:
    ch_bams  // channel: [ meta, path(bam), path(bai) ]
  main:
     ch_versions = Channel.empty()
     ch_germline_vcfs = Channel.empty()
      
    for ( tool in params.bam_germline_short_variant_discovery.tool ) {
      tool = tool.toLowerCase()      
      known_tool = false    

      if ( tool == "gatk4haplotypecaller" ) {
    
        def interval_list = file( params.genomes[params.genome].interval_list, checkIfExists: true )
        def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
        def fai = file( fasta.toString()+".fai", checkIfExists: true )
        def dict = file( fasta.toString().replace(".fasta",".dict"), checkIfExists: true )

        ch_intervals = Channel.value( interval_list )
          .map{ genome_interval -> [ [ id:'intervals' ], genome_interval ] }    
        ch_fasta = Channel.value( fasta )
          .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
        ch_fai = Channel.value( fai )
          .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }
        ch_dict = Channel.value( dict )
          .map{ genome_dict -> [ [ id:'dict' ], genome_dict ] }

        GATK4_SPLITINTERVALS( ch_intervals, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( GATK4_SPLITINTERVALS.out.versions )  

        ch_split_intervals = GATK4_SPLITINTERVALS.out.split_intervals
          .map{ meta, intervals -> [ intervals ]}
          .flatten()
                
        ch_bam_interval = ch_bams
          .combine( ch_split_intervals )
          .map{ meta, bam, bai, interval_file -> 
              m = interval_file =~ /(\d+)-scattered.interval_list/
              def interval = m[0][1]
              [ meta + [id: meta.id+"."+interval ], bam, bai, interval_file ]
            }

        BAM_BASE_RECALIBRATION_GATK4( ch_bam_interval, fasta, fai, dict )
        ch_versions = ch_versions.mix( BAM_BASE_RECALIBRATION_GATK4.out.versions ) 
        
        ch_bam_bqsr_interval = BAM_BASE_RECALIBRATION_GATK4.out.bam
          .combine( ch_split_intervals )
          .map{ meta, bam, bai, interval_file -> 
              m = interval_file =~ /(\d+)-scattered.interval_list/
              def interval = m[0][1]
              [ meta + [id: meta.id+"."+interval ], bam, bai, interval_file, [] ]
            }
            
        BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER( ch_bam_bqsr_interval, fasta, fai, dict )
        ch_versions = ch_versions.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.versions ) 

        ch_germline_vcfs = ch_germline_vcfs.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.vcf.join( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.tbi ) )
        
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_germline_vcfs // channel: [ meta, vcf, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


