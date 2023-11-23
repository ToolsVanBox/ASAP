//
// BAM GERMLINE SHORT VARIANT DISCOVERY HAPLOTYPECALLER

// Include nf-core modules
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_GENOTYPEGVCFS} from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_VQSR } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_COMBINEGVCFS } from '../../../modules/nf-core/gatk4/combinegvcfs/main'                                                
include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'                                              
include { GATK4_GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { BCFTOOLS_SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP } from '../../../modules/nf-core/gatk4/applyvqsr/main'

workflow BAM_GERMLINE_SHORT_VARIANT_DISCOVERY_GATK4HAPLOTYPECALLER {
  take:
    ch_bam_interval   // channel: [ val(meta), path(bam), path(bai), path(interval), path(dragstr_model) ]
    ch_fasta   // channel: [ val(meta), path(fasta) ]
    ch_fai   // channel: [ val(meta), path(fai) ]
    ch_dict   // channel: [ val(meta), path(dict) ]

  main:
    ch_versions = Channel.empty()
    ch_vcf = Channel.empty()
    ch_tbi = Channel.empty()
    
    fasta = ch_fasta.map{ meta, fasta -> [ fasta ] }
    fai = ch_fai.map{ meta, fai -> [ fai ] }
    dict = ch_dict.map{ meta, dict -> [ dict ] }
    
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

    // GATK4_HAPLOTYPECALLER.out.vcf.view()
    // IF ONE MAPPER AND MARKDUP
    ch_genomicsdbimport = GATK4_HAPLOTYPECALLER.out.vcf
      .join( GATK4_HAPLOTYPECALLER.out.tbi )
      .join( ch_bam_interval )
      .groupTuple( by: 5 )
      .map{ meta, gvcf, tbi, bam, bai, intervals, x ->
        m = meta[0].id =~ /.(\d+)$/
        def interval_id = m[0][1]
        [ [id: "${meta[0].run_id}.${interval_id}", run_id: meta[0].run_id ], gvcf, tbi, intervals, [], [] ]
      }
    //   GATK4_GENOMICSDBIMPORT( ch_genomicsdbimport, false, false, false )

    println "1"
    ch_genomicsdbimport.view()

      // IF MULTIPLE MAPPERS AND/OR MARKDUP
    ch_genomicsdbimport = ch_genomicsdbimport
      .map{ meta, vcf, tbi, intervals, x, y ->
        vcfs = vcf.toList()
        m = meta.id =~ /.(\d+)$/
        def interval_id = m[0][1]
        for ( v in vcfs ) {
          sample_map_file = file("sample_map.${interval_id}.txt")
          // if ( ! sample_map_file.exists() ) {
            sample_map_file.append( v.getName().toString().replaceAll(/.\d+.vcf.gz$/,"") + "\t" + v + "\n" )          
          // }
        }
        [ meta, sample_map_file, [], intervals, [], [] ]
      }

    println "2"
    ch_genomicsdbimport.view()

    GATK4_GENOMICSDBIMPORT( ch_genomicsdbimport, false, false, true )

    ch_versions = ch_versions.mix( GATK4_GENOMICSDBIMPORT.out.versions )  

    ch_genotypegvcfs = GATK4_GENOMICSDBIMPORT.out.genomicsdb
      .map{ meta, genomicsdb -> 
        [ meta, genomicsdb, [], [], [] ] 
      }
    
    GATK4_GENOTYPEGVCFS( ch_genotypegvcfs, fasta, fai, dict, dbsnp, dbsnp_tbi )
    ch_versions = ch_versions.mix( GATK4_GENOTYPEGVCFS.out.versions )

    BCFTOOLS_SORT( GATK4_GENOTYPEGVCFS.out.vcf )
    ch_versions = ch_versions.mix( BCFTOOLS_SORT.out.versions )


    ch_mergevcfs_genotypegvcfs = BCFTOOLS_SORT.out.vcf
      .map{ meta, vcf -> 
        [ [ id:meta.run_id, calling_type:"germline", caller:"gatk4haplotypecaller" ], vcf ]
      }
      .groupTuple()

    GATK4_MERGEVCFS_GENOTYPEGVCFS( ch_mergevcfs_genotypegvcfs, ch_dict )
    ch_versions = ch_versions.mix( GATK4_MERGEVCFS_GENOTYPEGVCFS.out.versions )
    
    known_indels = []
    known_indels_tbi = []
    known_indels_labels = []

    known_snps = []
    known_snps_tbi = []
    known_snps_labels = []

    for ( known_site_file in params.genomes[params.genome].known_sites ) {
      def known_site_tbi_file = known_site_file+".tbi"
      def file_exists = file(known_site_file).exists()
      def tbi_exists = file(known_site_tbi_file).exists()
      if (! (file_exists | tbi_exists) ) {
        println ">>> WARNING: Known sites file ${known_site_tbi_file} and/or ${known_site_tbi_file} does not exist and will not be included"
      }
      else {
        if ( known_site_file =~/indels/ ) {
          known_indels.add( known_site_file )
          known_indels_tbi.add( known_site_tbi_file )
          resource = file(known_site_file).getSimpleName()
          label = "--resource:${resource},${params.bam_germline_short_variant_discovery_haplotypecaller.gatk4variantrecalibrator_label} ${known_site_file}"
          known_indels_labels.add( label )
        }
        if ( known_site_file =~/snps/ ) {
          known_snps.add( known_site_file )
          known_snps_tbi.add( known_site_tbi_file )
          resource = file(known_site_file).getSimpleName()
          label = "--resource:${resource},${params.bam_germline_short_variant_discovery_haplotypecaller.gatk4variantrecalibrator_label} ${known_site_file}"
          known_snps_labels.add( label )
        }  
      }
    }

    ch_vqsr = GATK4_MERGEVCFS_GENOTYPEGVCFS.out.vcf
      .join( GATK4_MERGEVCFS_GENOTYPEGVCFS.out.tbi )
    
    GATK4_VARIANTRECALIBRATOR_INDEL( ch_vqsr, known_indels, known_indels_tbi, known_indels_labels, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_VARIANTRECALIBRATOR_INDEL.out.versions )

    GATK4_VARIANTRECALIBRATOR_SNP( ch_vqsr, known_snps, known_snps_tbi, known_snps_labels, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_VARIANTRECALIBRATOR_SNP.out.versions )


    ch_applyvqsr_snp = ch_vqsr
      .join( GATK4_VARIANTRECALIBRATOR_SNP.out.recal )
      .join( GATK4_VARIANTRECALIBRATOR_SNP.out.idx )
      .join( GATK4_VARIANTRECALIBRATOR_SNP.out.tranches )

    GATK4_APPLYVQSR_SNP( ch_applyvqsr_snp, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_APPLYVQSR_SNP.out.versions )

    ch_applyvqsr_indel = GATK4_APPLYVQSR_SNP.out.vcf
      .join( GATK4_APPLYVQSR_SNP.out.tbi )
      .join( GATK4_VARIANTRECALIBRATOR_INDEL.out.recal )
      .join( GATK4_VARIANTRECALIBRATOR_INDEL.out.idx )
      .join( GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches )
    
    GATK4_APPLYVQSR_INDEL( ch_applyvqsr_indel, fasta, fai, dict )
    ch_versions = ch_versions.mix( GATK4_APPLYVQSR_INDEL.out.versions )

    ch_merge_vqsr = GATK4_APPLYVQSR_INDEL.out.vcf
      .groupTuple()
    
    // ch_merge_vqsr = GATK4_APPLYVQSR_SNP.out.vcf
    //   .mix( GATK4_APPLYVQSR_INDEL.out.vcf )     
    //   .groupTuple()

    GATK4_MERGEVCFS_VQSR( ch_merge_vqsr, ch_dict )
    ch_versions = ch_versions.mix( GATK4_MERGEVCFS_VQSR.out.versions )  

    ch_vcf = ch_vcf.mix( GATK4_MERGEVCFS_VQSR.out.vcf )
    ch_tbi = ch_tbi.mix( GATK4_MERGEVCFS_VQSR.out.tbi )

  emit:
    vcf = ch_vcf // channel: [ meta, vcf ] 
    tbi = ch_tbi // channel: [ meta, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}


