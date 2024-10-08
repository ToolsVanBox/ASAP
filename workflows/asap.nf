
nextflow.enable.dsl=2

include { fromSamplesheet } from 'plugin/nf-validation'

// Include local subworkflows
include { FASTQ_QC_PRE_MAPPING } from '../subworkflows/local/fastq_qc_pre_mapping/main.nf'

include { FASTQ_ALIGN } from '../subworkflows/local/fastq_align/main.nf'
include { BAM_MERGE } from '../subworkflows/local/bam_merge/main.nf'
include { BAM_MARKDUPLICATES } from '../subworkflows/local/bam_markduplicates/main.nf'
include { BAM_QC_POST_MAPPING } from '../subworkflows/local/bam_qc_post_mapping/main.nf'
include { BAM_FINGERPRINT } from '../subworkflows/local/bam_fingerprint/main.nf'
include { BAM_CONVERT_TO_CRAM } from '../subworkflows/local/bam_convert_to_cram/main.nf'
include { BAM_TELOMERES } from '../subworkflows/local/bam_telomeres/main.nf'

include { BAM_GERMLINE_SHORT_VARIANT_DISCOVERY } from '../subworkflows/local/bam_germline_short_variant_discovery/main.nf'
include { BAM_GERMLINE_COPY_NUMBER_DISCOVERY } from '../subworkflows/local/bam_germline_copy_number_discovery/main.nf'
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_germline_structural_variant_discovery/main.nf'

include { BAM_SOMATIC_SHORT_VARIANT_DISCOVERY } from '../subworkflows/local/bam_somatic_short_variant_discovery/main.nf'
include { BAM_SOMATIC_COPY_NUMBER_DISCOVERY } from '../subworkflows/local/bam_somatic_copy_number_discovery/main.nf'
include { BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_somatic_structural_variant_discovery/main.nf'

include { BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_tumoronly_structural_variant_discovery/main.nf'
include { BAM_TUMORONLY_COPY_NUMBER_DISCOVERY } from '../subworkflows/local/bam_tumoronly_copy_number_discovery/main.nf'

include { VCF_GERMLINE_SHORT_VARIANT_FILTRATION } from '../subworkflows/local/vcf_germline_short_variant_filtration/main.nf'
include { VCF_SOMATIC_SHORT_VARIANT_FILTRATION } from '../subworkflows/local/vcf_somatic_short_variant_filtration/main.nf'
include { VCF_SHORT_VARIANT_ANNOTATION } from '../subworkflows/local/vcf_short_variant_annotation/main.nf'
include { VCF_STRUCTURAL_VARIANT_FILTRATION } from '../subworkflows/local/vcf_structural_variant_filtration/main.nf'

include { VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION } from '../subworkflows/local/vcf_germline_short_variant_somatic_filtration/main.nf'

include { BAM_HLA_TYPE_CALLING } from '../subworkflows/local/bam_hla_type_calling/main.nf'

 
// Include nf-core modules
include { GATK4_SPLITINTERVALS } from '../modules/nf-core/gatk4/splitintervals/main'                                                                          

ch_input = parseSampleSheet( Channel.fromSamplesheet("input") )

ch_input_type = ch_input.branch{
        bam:   it[0].data_type == "bam"
        fastq: it[0].data_type == "fastq"
        cram: it[0].data_type == "cram"
    }

workflow ASAP {

    // Create empty channels
    ch_versions = Channel.empty()
    ch_bam_bai = Channel.empty()
    ch_other_bam_bai = Channel.empty()
    ch_germline_vcfs = Channel.empty()
    ch_germline_tbi = Channel.empty()
    ch_somatic_vcfs = Channel.empty()
    ch_somatic_tbi = Channel.empty() 
    ch_input_annotation_vcfs = Channel.empty()
    ch_annotation_vcf_tbi = Channel.empty()
    ch_hla_vcf = Channel.empty()
    ch_hla_sc_vcf = Channel.empty()
    ch_hla_sf_vcf = Channel.empty()
    ch_rna_bam = Channel.empty()
    ch_cnv_dir = Channel.empty()
    
    // Define variables
    def fasta = file( params.genomes[params.genome].fasta, checkIfExists: true )
    def fai = file( params.genomes[params.genome].fai, checkIfExists: true )
    def dict = file( params.genomes[params.genome].dict, checkIfExists: true )
    def fasta_dict = file( params.genomes[params.genome].fasta_dict, checkIfExists: true )
    
    // Create channels of the variables
    ch_fasta = Channel.value( fasta )
      .map{ genome_fasta -> [ [ id:'fasta' ], genome_fasta ] }    
    ch_fai = Channel.value( fai )
      .map{ genome_fai -> [ [ id:'fai' ], genome_fai ] }
    ch_dict = Channel.value( dict )
      .map{ genome_dict -> [ [ id:'dict' ], genome_dict ] }
    ch_fasta_dict = Channel.value( dict )
      .map{ genome_fasta_dict -> [ [ id:'fasta_dict' ], genome_fasta_dict ] }

    // Split bams from the input data into dedupped and not dedupped (other)
    ch_bam_types = ch_input_type.bam.branch{
        dedup: it[1].toString() =~ /.*dedup.bam$/
        other: true
    }

    // Split crams from the input data into dedupped and not dedupped (other)
    ch_cram_types = ch_input_type.cram.branch{
        dedup: it[1].toString() =~ /.*dedup.cram$/
        other: true
    }

    // Assign the other type of bam and cram files 
    ch_other_bam_bai = ch_bam_types.other
    //ch_other_cram_crai = ch_cram_types.other

    // FASTQ QC PRE MAPPING
    if ( params.run.fastq_qc_pre_mapping ) {
        FASTQ_QC_PRE_MAPPING( ch_input_type.fastq )    
        ch_versions = ch_versions.mix( FASTQ_QC_PRE_MAPPING.out.versions )
    }

    // FASTQ ALIGN
    if ( params.run.fastq_align ) {
        FASTQ_ALIGN( ch_input_type.fastq, ch_fasta )
        ch_versions = ch_versions.mix( FASTQ_ALIGN.out.versions )

        // Create channel of bams to merge together
        ch_bams_to_merge = FASTQ_ALIGN.out.bam
            .map{ meta, bam ->
                meta = meta - meta.subMap('read_group')
                meta.id = meta.sample+"."+meta.align
                [ meta, bam ]
            }
            .groupTuple()

        // MERGE BAM FILES
        BAM_MERGE( ch_bams_to_merge, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_MERGE.out.versions.first() )

        // Add non dedupped bams files to the "other" channel
        ch_other_bam_bai = ch_other_bam_bai.mix( BAM_MERGE.out.bam.join( BAM_MERGE.out.bai ) )
    }
    
    // MARK DUPLICATES
    if ( params.run.bam_markduplicates ) {

        BAM_MARKDUPLICATES( ch_other_bam_bai, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_MARKDUPLICATES.out.versions )

        // Join bam files and bai files    
        ch_dedup_bam_bai = BAM_MARKDUPLICATES.out.bam
            .join( BAM_MARKDUPLICATES.out.bai )
        
        // Mix dedupped input bam files with the "new" dedupped bam files
        ch_bam_bai = ch_bam_types.dedup
            .mix( ch_dedup_bam_bai )
    } else {        
        // Mix dedupped input bam files with the non dedupped bam files IF Mark duplicated is skipped 
        ch_bam_bai = ch_bam_types.dedup
            .mix( ch_other_bam_bai )
    }

    // Convert bam files to cram files (for backup)
    if ( params.run.bam_convert_to_cram ) {
        BAM_CONVERT_TO_CRAM( ch_bam_bai, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_CONVERT_TO_CRAM.out.versions.first() )
        ch_cram_crai = BAM_CONVERT_TO_CRAM.out.cram_crai
    }

    // POST MAPPING
    if ( params.run.bam_qc_post_mapping ) {
        BAM_QC_POST_MAPPING( ch_bam_bai, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_QC_POST_MAPPING.out.versions )
    }

    // FINGERPRINT
    if ( params.run.bam_fingerprint ) {
        BAM_FINGERPRINT( ch_bam_bai, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_FINGERPRINT.out.versions )
    }
    
    // GET TELOMERES
    if ( params.run.bam_telomeres ) {
        BAM_TELOMERES( ch_bam_bai )
        ch_versions = ch_versions.mix( BAM_TELOMERES.out.versions )
    }   
   
    // Variant calling
    if ( params.run.bam_germline_short_variant_discovery || params.run.bam_somatic_short_variant_discovery ) {
        def interval_list = file( params.genomes[params.genome].interval_list, checkIfExists: true )
    
        ch_intervals = Channel.value( interval_list )
            .map{ genome_interval -> [ [ id:'intervals' ], genome_interval ] }    

        GATK4_SPLITINTERVALS( ch_intervals, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( GATK4_SPLITINTERVALS.out.versions )  

        ch_split_intervals = GATK4_SPLITINTERVALS.out.split_intervals
            .map{ meta, intervals -> 
                [ intervals ]
            }
            .flatten()
            .map{ interval_file ->
                m = file(interval_file).getName() =~ /^(\d+)-scattered.interval_list/
                [ [id: m[0][1] ], interval_file]
            }
    }

    // Germline short variant calling
    if ( params.run.bam_germline_short_variant_discovery ) {        
        BAM_GERMLINE_SHORT_VARIANT_DISCOVERY( ch_bam_bai, ch_split_intervals, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY.out.versions )

        ch_germline_vcfs = BAM_GERMLINE_SHORT_VARIANT_DISCOVERY.out.vcf 
        ch_germline_tbi = BAM_GERMLINE_SHORT_VARIANT_DISCOVERY.out.tbi 
        // Variant filtration
        if ( params.run.vcf_germline_short_variant_filtration ) {
            ch_germline_variant_filtration = ch_germline_vcfs
                .join( ch_germline_tbi )

            VCF_GERMLINE_SHORT_VARIANT_FILTRATION( ch_germline_variant_filtration, ch_fasta, ch_fai, ch_dict )
            ch_versions = ch_versions.mix( VCF_GERMLINE_SHORT_VARIANT_FILTRATION.out.versions )

            ch_germline_vcfs = VCF_GERMLINE_SHORT_VARIANT_FILTRATION.out.vcf
            ch_germline_tbi = VCF_GERMLINE_SHORT_VARIANT_FILTRATION.out.tbi
        }
    }

    // Germline copy number discovery 
    if ( params.run.bam_germline_copy_number_discovery ) {
        BAM_GERMLINE_COPY_NUMBER_DISCOVERY( ch_bam_bai, ch_fasta )
        ch_versions = ch_versions.mix( BAM_GERMLINE_COPY_NUMBER_DISCOVERY.out.versions )
    }   

    // Germline structural variant discovery 
    if ( params.run.bam_germline_structural_variant_discovery ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY( ch_bam_bai, ch_fasta, ch_fai, ch_fasta_dict )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY.out.versions )

        ch_gridss_vcf = BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY.out.gridss_vcf
    }

    ch_bam_bai_sample_type = ch_bam_bai.branch{
        normal: it[0].sample_type == "normal"
        tumor: it[0].sample_type == "tumor"
    }
     
    // Somatic calling input preparation 
    if ( params.run.bam_somatic_structural_variant_discovery || params.run.bam_somatic_copy_number_discovery || params.run.bam_somatic_short_variant_discovery ) {
        ch_bam_bai_normal_tumor = ch_bam_bai_sample_type.normal.ifEmpty( [ [], null, null ] )
                .combine( ch_bam_bai_sample_type.tumor.ifEmpty( [ [], null, null ] ) )
                .map{ normal_meta, normal_bam, normal_bai, tumor_meta, tumor_bam, tumor_bai ->
                    if ( normal_bam == null ) { error("Need a normal sample for somatic calling") }
                    if ( tumor_bam == null ) { error("Need a tumor sample for somatic calling") }
                    [ [ id: "${normal_meta.id}_${tumor_meta.id}", run_id: normal_meta.run_id, tumor_sample_id: tumor_meta.id, normal_sample_id: normal_meta.id ], normal_bam, normal_bai, tumor_bam, tumor_bai ]
                }       
    }

    // Somatic structural variant discovery 
    if ( params.run.bam_somatic_structural_variant_discovery ) {                    
        BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY( ch_bam_bai_normal_tumor, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY.out.versions )
    }

    // Somatic copy number discovery 
    if ( params.run.bam_somatic_copy_number_discovery ) {                    
        BAM_SOMATIC_COPY_NUMBER_DISCOVERY( ch_bam_bai_sample_type.normal, ch_bam_bai_sample_type.tumor, ch_fasta )
        ch_versions = ch_versions.mix( BAM_SOMATIC_COPY_NUMBER_DISCOVERY.out.versions )
    }

    // Somatic short variant discovery 
    if ( params.run.bam_somatic_short_variant_discovery ) {
        BAM_SOMATIC_SHORT_VARIANT_DISCOVERY( ch_bam_bai_sample_type.tumor, ch_bam_bai_sample_type.normal, ch_split_intervals, ch_fasta, ch_fai, ch_dict )
        ch_versions = ch_versions.mix( BAM_SOMATIC_SHORT_VARIANT_DISCOVERY.out.versions )
        
        ch_somatic_vcfs = BAM_SOMATIC_SHORT_VARIANT_DISCOVERY.out.vcf 
        ch_somatic_tbi = BAM_SOMATIC_SHORT_VARIANT_DISCOVERY.out.tbi 
        ch_somatic_f1r2 = BAM_SOMATIC_SHORT_VARIANT_DISCOVERY.out.f1r2 
        ch_somatic_stats = BAM_SOMATIC_SHORT_VARIANT_DISCOVERY.out.stats 

        // Variant filtration 
        if ( params.run.vcf_somatic_short_variant_filtration ) {
            VCF_SOMATIC_SHORT_VARIANT_FILTRATION(ch_somatic_vcfs, ch_somatic_tbi, ch_somatic_f1r2, ch_somatic_stats, ch_bam_bai_sample_type.tumor, ch_bam_bai_sample_type.normal, ch_split_intervals, ch_fasta, ch_fai, ch_dict  )
            ch_versions = ch_versions.mix( VCF_SOMATIC_SHORT_VARIANT_FILTRATION.out.versions )

            ch_somatic_vcfs = VCF_SOMATIC_SHORT_VARIANT_FILTRATION.out.vcf
            ch_somatic_tbi = VCF_SOMATIC_SHORT_VARIANT_FILTRATION.out.tbi
        }

        // Assign hla vcf to be the somatic vcf (Double check if this takes along the filtration step)
        ch_hla_sc_vcf = ch_somatic_vcfs.join(ch_somatic_tbi)
    }
    
    // Tumor-only
    if ( params.run.bam_tumoronly_structural_variant_discovery || params.run.bam_tumoronly_copy_number_discovery ) {    
        ch_bam_bai_tumor = ch_bam_bai_sample_type.tumor.ifEmpty( [ [], null, null ] )
            .map{ meta, bam, bai ->
                if ( bam == null ) { error("Need a tumor sample for tumor only calling") }
                [ meta, bam, bai ]
            }
    }

    // Tumor only structural variant discovery 
    if ( params.run.bam_tumoronly_structural_variant_discovery ) {                    
        BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY( ch_bam_bai_tumor, ch_fasta, ch_fai )
        ch_versions = ch_versions.mix( BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY.out.versions )
    }

    // Tumor only copy number discovery 
    if ( params.run.bam_tumoronly_copy_number_discovery ) {                    
        BAM_TUMORONLY_COPY_NUMBER_DISCOVERY( ch_bam_bai_tumor, ch_fasta )
        ch_versions = ch_versions.mix( BAM_TUMORONLY_COPY_NUMBER_DISCOVERY.out.versions )
    }
    
    // Short variant annotation 
    if ( params.run.vcf_short_variant_annotation ) {
        if (params.run.bam_germline_short_variant_discovery){
            ch_input_annotation_vcfs = ch_input_annotation_vcfs.mix(ch_germline_vcfs)
        }
        if (params.run.bam_somatic_short_variant_discovery){
            ch_input_annotation_vcfs = ch_input_annotation_vcfs.mix(ch_somatic_vcfs)
        }

        VCF_SHORT_VARIANT_ANNOTATION( ch_input_annotation_vcfs, ch_fasta )
        ch_annotation_vcf_tbi = VCF_SHORT_VARIANT_ANNOTATION.out.vcf_tbi

        // Assign the hla vcf if you have a somatic output vcf
        ch_hla_sc_vcf = ch_annotation_vcf_tbi.map{ meta, vcf, tbi ->
            if ( meta.calling_type == "somatic" ) {
                [ meta, vcf, tbi]
            }
        }
    }

    // Structural variant filtration
    if ( params.run.vcf_somatic_structural_variant_filtration || params.run.vcf_germline_structural_variant_filtration || params.run.vcf_tumoronly_structural_variant_filtration ) {
        
        VCF_STRUCTURAL_VARIANT_FILTRATION( ch_gridss_vcf, ch_bam_bai_sample_type.normal, ch_bam_bai_sample_type.tumor, ch_fasta, ch_fai, ch_fasta_dict )
        ch_versions = ch_versions.mix( VCF_STRUCTURAL_VARIANT_FILTRATION.out.versions )
    }

    // Somatic filtering
    if ( params.run.vcf_germline_short_variant_somatic_filtration ) {
        ch_germline_vcf_tbi = ch_annotation_vcf_tbi
            .map{ meta, vcf, tbi ->
                if ( meta.calling_type == "germline" ) {
                    [ meta, vcf, tbi]
                }
            }
        
        VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION( ch_germline_vcf_tbi, ch_bam_bai )
        ch_versions = ch_versions.mix( VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION.out.versions )
        
        VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION.out.view()
        //ch_hla_sf_vcf = VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION.out.vcf.join( VCF_GERMLINE_SHORT_VARIANT_SOMATIC_FILTRATION.out.tbi )
    }

    // HLA type calling 
    if (params.run.bam_hla_type_calling) {
        ch_hla_vcf = ch_hla_sc_vcf.mix( ch_hla_sf_vcf )

        BAM_HLA_TYPE_CALLING(ch_bam_bai_sample_type.normal, ch_bam_bai_sample_type.tumor, ch_rna_bam, ch_cnv_dir, ch_hla_vcf,  ch_fasta, ch_fai)
        ch_versions = ch_versions.mix( BAM_HLA_TYPE_CALLING.out.versions )
    }

    
//     // CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
//     // version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
}

def parseSampleSheet( ch_csv ) {
    if ( ! params.run_id ) {
        run_id = params.out_dir.split('/')[-1]
    if ( params.out_dir == './' ) {
        run_id = "${PWD}".split('/')[-1]
    }
    } else {
        run_id = params.run_id
    }

    ch_csv
        .map{ meta, fastq_1, fastq_2, bam, bai, sample_type, platform ->
            meta = meta + [ run_id: run_id, sample_type: sample_type ]
            if ( fastq_1 && fastq_2 ) {
                // Map per lane OLD WAY

                // def (flowcell, lane, machine, run_nr) = flowcellLaneFromFastq(fastq_1)
                // def rg_id = "${meta.sample}_${flowcell}_${lane}"
                // def read_group  = "\"@RG\\tID:${rg_id}\\tSM:${meta.sample}\\tPL:${platform}\\tLB:${meta.sample}\\tPU:${flowcell}\""
                
                // Map per fastq NEW APPROACH                
                def rg_id = fastq_1.getSimpleName()
                def read_group  = "\"@RG\\tID:${rg_id}\\tSM:${meta.sample}\\tPL:${platform}\\tLB:${meta.sample}\""
                meta = meta + [ id: rg_id.toString() ]
                meta = meta + [ read_group: read_group.toString() ]
                meta = meta + [ data_type: "fastq" ]
                return [ meta, [ fastq_1, fastq_2 ] ]
            }
            if ( bam ) {
                meta = meta + [ id: bam.getBaseName() ]
                if ( bam.getExtension() == "bam" ) {
                    meta = meta + [ data_type: "bam" ]
                }
                if ( bam.getExtension() == "cram" ) {
                    meta = meta + [ data_type: "cram" ]
                }
                return [ meta, bam, bai ]
            }
            
        }
}

def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        machine = fields[0]
        run_nr = fields[1].toInteger()
        fcid = fields[2]
        lane = fields[3].toInteger()

    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    return [fcid, lane, machine, run_nr]
}
