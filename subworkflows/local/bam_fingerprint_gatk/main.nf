//
// BAM FINGERPRINT GATK

// Include nf-core modules
include { GATK_UNIFIEDGENOTYPER } from '../../../modules/nf-core/gatk/unifiedgenotyper/main'
include { TABIX_BGZIP } from '../../../modules/nf-core/tabix/bgzip/main'
include { FINGERPRINT_HEATMAP } from '../../../modules/local/fingerprint/heatmap/main'

workflow BAM_FINGERPRINT_GATK {
	take: 
		ch_bams		// channel: [ val(meta), path(bam), path(bai) ]
		ch_fasta	// channel: [ val(meta), path(fasta) ]
		ch_fai		// channel: [ val(meta), path(fai) ]
		ch_dict		// channel: [ val(meta), path(dict) ]

	main: 
		ch_versions = Channel.empty()
		ch_fingerprintvcf = channel.empty()

		def snpfile = file( params.genomes[params.genome].gatk_fingerprint_vcf, checkIfExists: true )
		ch_snpfile = Channel.value( snpfile )
	  		.map{ gatk_snpfile -> [ [ id:'snpfile' ], gatk_snpfile ] }   

	  	// Create vcf 
		GATK_UNIFIEDGENOTYPER( ch_bams, ch_fasta, ch_fai, ch_dict, ch_snpfile, [[],[]], [[],[]], [[],[]])
		ch_versions = ch_versions.mix( GATK_UNIFIEDGENOTYPER.out.versions )  

		// Unzip the vcf file 
		TABIX_BGZIP( GATK_UNIFIEDGENOTYPER.out.vcf.map{ meta, vcfs ->
				vcf = vcfs.find { it.toString().endsWith("gatk_unifiedgenotyper.vcf.gz")}
				[ meta, vcf]
			})
		ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )  

		// Combine vcf file
		ch_fingerprintvcf = TABIX_BGZIP.out.output
			.map{ meta, vcf ->
				[ [ id:meta.run_id ], vcf]
			}.groupTuple()

		// Add heatmap (check bam_copynumber freec) 
		FINGERPRINT_HEATMAP(ch_fingerprintvcf)

	emit:
		versions = ch_versions					   // channel: [ versions.yml ]
}







