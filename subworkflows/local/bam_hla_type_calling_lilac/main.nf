/*
 LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
*/

//import Constants
//import Utils

// Include local subworkflows
//include { CUSTOM_EXTRACTCONTIG as EXTRACTCONTIG } from '../../../modules/local/hmftools/lilac/lilac_extract_and_index_contig/main'
//include { CUSTOM_REALIGNREADS as REALIGNREADS   } from '../../../modules/local/hmftools/lilac/lilac_realign_reads_lilac/main'
//include { CUSTOM_SLICE as SLICEBAM              } from '../../../modules/local/hmftools/lilac/lilac_slice/main'
include { LILAC                                 } from '../../../modules/local/hmftools/lilac/lilac_calling/main'


workflow BAM_HLA_TYPE_CALLING_LILAC {
    take: 
        // Sample data
        ch_input           // channel: [mandatory] [ meta ]
        ch_normal_bam      // channel: [mandatory] [ meta, bam, bai ]
        ch_tumor_bam       // channel: [mandatory] [ meta, bam, bai ]
        ch_tumor_rna_bam   // channel: [mandatory] [ meta, bam, bai ]
        ch_cnv_dir         // channel: [mandatory] [ meta, cnv_dir ]
        ch_somatic_vcf     // channel: [mandatory] [ meta, path(vcf) ]

        // Reference data
        genome_fasta       // channel: [mandatory] /path/to/genome_fasta
        genome_fai         // channel: [mandatory] /path/to/genome_fai


    main: 
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Parse the fixed variable to be taken by lilac 
        genome_fasta = genome_fasta.map{meta, fasta_path -> [ fasta_path ] }.collect()

        // Do the optionial parameter check 
        normal_dna_bam_ch = ch_normal_bam ? ch_normal_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty() 
        normal_dna_bai_ch = ch_normal_bam ? ch_normal_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        tumor_dna_bam_ch = ch_tumor_bam ? ch_tumor_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty()
        tumor_dna_bai_ch = ch_tumor_bam ? ch_tumor_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        tumor_rna_bam_ch = ch_tumor_rna_bam ? ch_tumor_rna_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty()
        tumor_rna_bai_ch = ch_tumor_rna_bam ? ch_tumor_rna_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        purple_dir_ch = ch_cnv_dir ? ch_cnv_dir.map{meta, cnv_dir -> [ cnv_dir ] } : Channel.empty()
        somatic_vcf_ch = ch_somatic_vcf ? ch_somatic_vcf.map{meta, vcf, vcf_tbi -> [ vcf ] } : Channel.empty()

        // Run Lilac with optional parameters 
        LILAC(
            ch_input, 
            normal_dna_bam_ch.ifEmpty([]), 
            normal_dna_bai_ch.ifEmpty([]),
            tumor_dna_bam_ch.ifEmpty([]),
            tumor_dna_bai_ch.ifEmpty([]),
            tumor_rna_bam_ch.ifEmpty([]), 
            tumor_rna_bai_ch.ifEmpty([]), 
            purple_dir_ch.ifEmpty([]), 
            somatic_vcf_ch.ifEmpty([]),
            genome_fasta,
            params.genomes[params.genome].lilac_chr_version, 
            params.genomes[params.genome].lilac_dir
        )

        // Add the info to the versions run in the pipeline 
        ch_versions = ch_versions.mix(LILAC.out.versions)

    emit: 
        lilac_dir = LILAC.out.lilac_dir  // channel: [ meta, lilac_dir ]

        versions  = ch_versions // channel: [ versions.yml ]
}