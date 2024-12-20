/*
 LILAC is a WGS tool for HLA typing and somatic CNV and SNV calling
*/

//import Constants
//import Utils

// Include local subworkflows
include { CUSTOM_EXTRACTCONTIG as EXTRACTCONTIG } from '../../../modules/local/hmftools/lilac/lilac_extract_and_index_contig/main'
include { CUSTOM_REALIGNREADS as REALIGNREADS   } from '../../../modules/local/hmftools/lilac/lilac_realign_reads_lilac/main'
include { CUSTOM_SLICE as SLICEBAM              } from '../../../modules/local/hmftools/lilac/lilac_slice/main'
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
        hla_slice_bed      // channel: [mandatory] /path/to/hla_slice_bed


    main: 
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()
        ch_slice_input = Channel.empty()
        ch_realign_inputs_sorted = Channel.empty()

        // Parse the fixed variable to be taken by lilac 
        genome_fasta = genome_fasta.map{meta, fasta_path -> [ fasta_path ] }.collect()
        genome_fai = genome_fai.map{meta, fai_path -> [ fai_path ] }.collect()

        // Do the optionial parameter check 
        tumor_rna_bam_ch = ch_tumor_rna_bam ? ch_tumor_rna_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty()
        tumor_rna_bai_ch = ch_tumor_rna_bam ? ch_tumor_rna_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        purple_dir_ch = ch_cnv_dir ? ch_cnv_dir.map{meta, cnv_dir -> [ cnv_dir ] } : Channel.empty()
        somatic_vcf_ch = ch_somatic_vcf ? ch_somatic_vcf.map{meta, vcf, vcf_tbi -> [ vcf ] } : Channel.empty()

        // Add the alt genome part here 
        if ( params.genomes[params.genome].lilac_genome_type == "alt" ){
            // Parse input channel
            ch_realign_inputs_sorted = ch_realign_inputs_sorted.mix(ch_tumor_bam, ch_normal_bam)

            // Custom bam slice 
            SLICEBAM(
                ch_realign_inputs_sorted,
                hla_slice_bed,
            )
            ch_versions = ch_versions.mix(SLICEBAM.out.versions)
            
            // Extract contig 
            ch_extract_contig_run = ch_realign_inputs_sorted
                .toList()
                .map { !it.isEmpty() }

            EXTRACTCONTIG(
                params.genomes[params.genome].lilac_chr,
                genome_fasta,
                genome_fai,
                ch_extract_contig_run,
            )
            ch_versions = ch_versions.mix(EXTRACTCONTIG.out.versions)
            
            // Realign reads
            ch_realign_bam = SLICEBAM.out.bam.map{ meta, bam, bai -> 
                meta.sample_id = meta.id
                [meta, bam, bai ]
            }
            
            REALIGNREADS(
                ch_realign_bam,
                EXTRACTCONTIG.out.contig,
                EXTRACTCONTIG.out.bwamem2_index,
            )
            ch_versions = ch_versions.mix(REALIGNREADS.out.versions)

            // Put to channels
            ch_slice_reunited_bams = REALIGNREADS.out.bam.branch { meta, bam, bai ->
                tumor: meta.sample_type == 'tumor'
                normal: meta.sample_type == 'normal'
            }

            // Fill variables 
            normal_dna_bam_ch = ch_slice_reunited_bams.normal ? ch_slice_reunited_bams.normal.map{ meta, bam, bai -> [ bam ] } : Channel.empty() 
            normal_dna_bai_ch = ch_slice_reunited_bams.normal ? ch_slice_reunited_bams.normal.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
            tumor_dna_bam_ch = ch_slice_reunited_bams.tumor ? ch_slice_reunited_bams.tumor.map{ meta, bam, bai -> [ bam ] } : Channel.empty() 
            tumor_dna_bai_ch = ch_slice_reunited_bams.tumor ? ch_slice_reunited_bams.tumor.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        } else {
            // Without alt genomes, use the default genomes 
            normal_dna_bam_ch = ch_normal_bam ? ch_normal_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty() 
            normal_dna_bai_ch = ch_normal_bam ? ch_normal_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
            tumor_dna_bam_ch = ch_tumor_bam ? ch_tumor_bam.map{ meta, bam, bai -> [ bam ] } : Channel.empty()
            tumor_dna_bai_ch = ch_tumor_bam ? ch_tumor_bam.map{ meta, bam, bai -> [ bai ] } : Channel.empty()
        }
        

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