//
// BAM GERMLINE SHORT VARIANT DISCOVERY HAPLOTYPECALLER

// Include nf-core modules
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS } from '../../../modules/nf-core/gatk4/mergevcfs/main'                                                      
include { GATK4_COMBINEGVCFS } from '../../../modules/nf-core/gatk4/combinegvcfs/main'                                                
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'                                              


workflow BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER {
  take:
    ch_bam_interval   // channel: [ val(meta), path(bam), path(bai), path(interval), path(dragstr_model) ]
    fasta   // path(fasta)
    fai   // path(fai)
    dict   // path(dict)
  main:
    ch_versions = Channel.empty()
    
    ch_dict = Channel.value( dict )
      .map{ genome_dict -> [ [ id:'dict' ], genome_dict ] }  
    
    dbsnp = []
    dbsnp_tbi = []

    for ( dbsnp_file in params.genomes[params.genome].dbsnp ) {
      def dbsnp_tbi_file = dbsnp_file+".tbi"
      def file_exists = file(dbsnp_file).exists()
      def tbi_exists = file(dbsnp_tbi_file).exists()
      if (! (file_exists | tbi_exists) ) {
        println ">>> WARNING: DBSNP file ${dbsnp_file} and/or ${dbsnp_tbi_file} does not exist and will not be included"
      }
      else {
        dbsnp.add( dbsnp_file )
        dbsnp_tbi.add( dbsnp_tbi_file )
      }
    }

    GATK4_HAPLOTYPECALLER( ch_bam_interval, fasta, fai, dict, dbsnp, dbsnp_tbi )
    ch_versions = ch_versions.mix( GATK4_HAPLOTYPECALLER.out.versions )  

    ch_gvcfs = GATK4_HAPLOTYPECALLER.out.vcf
      .map{ meta, vcf ->
        [ [id: meta.sample_id, run_id: meta.run_id], vcf ]
      }
      .groupTuple()

    GATK4_MERGEVCFS( ch_gvcfs, ch_dict )
    ch_versions = ch_versions.mix( GATK4_MERGEVCFS.out.versions )  

    ch_gvcfs = GATK4_MERGEVCFS.out.vcf
      .join( GATK4_MERGEVCFS.out.tbi)
      .map{ meta, vcf, tbi ->
        [ [id:meta.run_id], vcf, tbi]
      }
      .groupTuple()
    
    GATK4_COMBINEGVCFS( ch_gvcfs, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_COMBINEGVCFS.out.versions )  

    ch_combined_gvcfs = GATK4_COMBINEGVCFS.out.combined_gvcf
      .join( GATK4_COMBINEGVCFS.out.combined_tbi)
      .map{ meta, gvcf, tbi ->
        [ meta, gvcf, tbi, [], [] ]
      }

    GATK4_GENOTYPEGVCFS( ch_combined_gvcfs, fasta, fai, dict, dbsnp, dbsnp_tbi )
    ch_versions = ch_versions.mix( GATK4_GENOTYPEGVCFS.out.versions )  

    GATK4_GENOTYPEGVCFS.out.vcf.view()
    GATK4_GENOTYPEGVCFS.out.tbi.view()

  emit:
    vcf = GATK4_GENOTYPEGVCFS.out.vcf // channel: [ meta, vcf ] 
    tbi = GATK4_GENOTYPEGVCFS.out.tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


