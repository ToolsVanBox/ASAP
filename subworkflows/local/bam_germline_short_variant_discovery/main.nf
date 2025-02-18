//
// BAM GERMLINE SHORT VARIANT DISCOVERY

// Include local subworkflows
include { BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER } from '../../../subworkflows/local/bam_germline_short_variant_discovery_gatk4haplotypecaller/main'       
include { BAM_BASE_RECALIBRATION_GATK4 } from '../../../subworkflows/local/bam_base_recalibration_gatk4/main'

// Include nf-core modules
include { GATK4_SPLITINTERVALS } from '../../../modules/nf-core/gatk4/splitintervals/main'                                                                          

workflow BAM_GERMLINE_SHORT_VARIANT_DISCOVERY {
  take:
    ch_bam_bai  // channel: [ meta, path(bam), path(bai) ]
    ch_intervals // channel: [ val(meta), path(interval_list) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]
    ch_dict // channel: [ val(meta), path(dict) ]

  main:
     ch_versions = Channel.empty()
     ch_germline_vcfs = Channel.empty()
     ch_germline_tbi = Channel.empty()
      
    for ( tool in params.bam_germline_short_variant_discovery.tool ) {
      def tool = tool.toLowerCase()      
      def known_tool = false    

      if ( tool == "gatk4haplotypecaller" ) {
                       
        def ch_bam_bai_interval = ch_bam_bai
          .combine( ch_intervals )
          .map{ meta, bam, bai, meta2, interval_file -> 
              meta = meta + [ calling_type: "germline" ] 
              meta = meta + [ short_variant_caller: "gatk4haplotypecaller" ]
              meta = meta + [ interval: meta2.id ] 
              meta.id = meta.id+".germline.gatk4haplotypecaller."+meta2.id
              [ meta, bam, bai, interval_file ]
            }

        BAM_BASE_RECALIBRATION_GATK4( ch_bam_bai_interval, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_BASE_RECALIBRATION_GATK4.out.versions ) 
        
        def ch_bam_bai_interval_bqsr = BAM_BASE_RECALIBRATION_GATK4.out.bam
          .join( BAM_BASE_RECALIBRATION_GATK4.out.bai )
          .combine( ch_intervals )
          .map{ meta, bam, bai, meta2, interval_file ->
              meta = meta + [ interval: meta2.id ] 
              meta.id = meta.id+".germline.gatk4haplotypecaller."+meta2.id
              [ meta, bam, bai, interval_file, [] ]
            }

        BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER( ch_bam_bai_interval_bqsr, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.versions ) 

        ch_germline_vcfs = ch_germline_vcfs.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.vcf )
        ch_germline_tbi = ch_germline_tbi.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER.out.tbi )
        
        known_tool = true
      }

      if ( ! known_tool ) {
       println ("WARNING: Skip ${tool}, because it's not known as a bam germline short variant discovery tool, this tool is not build in (yet).")
      }
    }
  emit:
    vcf = ch_germline_vcfs // channel: [ meta, vcf ]
    tbi = ch_germline_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


