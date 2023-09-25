
nextflow.enable.dsl=2

// Include local subworkflows
include { FASTQ_QC_PRE_MAPPING } from '../subworkflows/local/fastq_qc_pre_mapping/main.nf'
include { FASTQ_ALIGN } from '../subworkflows/local/fastq_align/main.nf'
include { BAM_MARKDUPLICATES } from '../subworkflows/local/bam_markduplicates/main.nf'
include { BAM_QC_POST_MAPPING } from '../subworkflows/local/bam_qc_post_mapping/main.nf'
include { BAM_FINGERPRINT } from '../subworkflows/local/bam_fingerprint/main.nf'
include { BAM_GERMLINE_SHORT_VARIANT_DISCOVERY } from '../subworkflows/local/bam_germline_short_variant_discovery/main.nf'
include { BAM_GERMLINE_COPY_NUMBER_DISCOVERY } from '../subworkflows/local/bam_germline_copy_number_discovery/main.nf'
include { BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_germline_structural_variant_discovery/main.nf'
include { BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_somatic_structural_variant_discovery/main.nf'
include { BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY } from '../subworkflows/local/bam_tumoronly_structural_variant_discovery/main.nf'

// Include nf-core modules
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'

input_sample = parseSampleSheet( params.input )
input_sample_type = input_sample.branch{
        bam:   it[0].data_type == "bam"
        fastq: it[0].data_type == "fastq"
    }

workflow WAP {

    ch_versions = Channel.empty()
    ch_bams = Channel.empty()
    ch_other_bams = Channel.empty()
    ch_dedup_bams = Channel.empty()

    ch_bam_types = input_sample_type.bam.branch{
        dedup: it[1].toString() =~ /.*dedup.bam$/
        other: true
    }

    if ( params.run.fastq_qc_pre_mapping ) {
        FASTQ_QC_PRE_MAPPING( input_sample_type.fastq )    
        ch_versions = ch_versions.mix( FASTQ_QC_PRE_MAPPING.out.versions )
    }

    if ( params.run.fastq_align ) {
        FASTQ_ALIGN( input_sample_type.fastq )
        ch_versions = ch_versions.mix( FASTQ_ALIGN.out.versions )

        ch_other_bams = ch_bam_types.other.mix( FASTQ_ALIGN.out.bam )
    }
    
    if ( params.run.bam_markduplicates ) {
        BAM_MARKDUPLICATES( ch_other_bams )
        ch_versions = ch_versions.mix( BAM_MARKDUPLICATES.out.versions )

        ch_bams = ch_bam_types.dedup
            .mix( BAM_MARKDUPLICATES.out.bam)
            .map{ meta, bam ,bai -> 
                [[id: meta.id, sample_id: meta.sample_id, run_id: run_id, sample_type: meta.sample_type], bam ,bai ]
            }
    } else {
        ch_bams = ch_bam_types.dedup
            .mix( ch_other_bams )
            .map{ meta, bam ,bai -> 
                [[id: meta.id, sample_id: meta.sample_id, run_id: run_id, sample_type: meta.sample_type], bam ,bai ]
            }
    }

    ch_bams_sample_type = ch_bams.branch{
        normal: it[0].sample_type == "normal"
        tumor: it[0].sample_type == "tumor"
    }

    ch_bams_tumor_normal = ch_bams_sample_type.tumor
            .combine( ch_bams_sample_type.normal.ifEmpty( [ [], null, null ] ) )
    
    if ( params.run.bam_qc_post_mapping ) {
        BAM_QC_POST_MAPPING( ch_bams )
        ch_versions = ch_versions.mix( BAM_QC_POST_MAPPING.out.versions )
    }

    if ( params.run.bam_fingerprint ) {
        BAM_FINGERPRINT( ch_bams )
        ch_versions = ch_versions.mix( BAM_FINGERPRINT.out.versions )
    }
    
    if ( params.run.bam_germline_copy_number_discovery ) {
        BAM_GERMLINE_COPY_NUMBER_DISCOVERY( ch_bams )
        ch_versions = ch_versions.mix( BAM_GERMLINE_COPY_NUMBER_DISCOVERY.out.versions )
    }

    if ( params.run.bam_germline_short_variant_discovery ) {
        BAM_GERMLINE_SHORT_VARIANT_DISCOVERY( ch_bams )
        ch_versions = ch_versions.mix( BAM_GERMLINE_SHORT_VARIANT_DISCOVERY.out.versions )
    }

    if ( params.run.bam_germline_structural_variant_discovery ) {
        BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY( ch_bams )
        ch_versions = ch_versions.mix( BAM_GERMLINE_STRUCTURAL_VARIANT_DISCOVERY.out.versions )
    }

    if ( params.run.bam_somatic_structural_variant_discovery ) {                    
        BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY( ch_bams )
        ch_versions = ch_versions.mix( BAM_SOMATIC_STRUCTURAL_VARIANT_DISCOVERY.out.versions )
    }

    if ( params.run.bam_tumoronly_structural_variant_discovery ) {                    
        BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY( ch_bams )
        ch_versions = ch_versions.mix( BAM_TUMORONLY_STRUCTURAL_VARIANT_DISCOVERY.out.versions )
    }

    
    // CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    // version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
}

def parseSampleSheet(csv_file) {

    if ( ! params.run_id ) {
    run_id = params.out_dir.split('/')[-1]
    if ( params.out_dir == './' ) {
        run_id = "${PWD}".split('/')[-1]
    }
    } else {
        run_id = params.run_id
    }

    Channel.fromPath(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.sample_id)) {
                error("Missing field in csv file header. The csv file must have field 'sample_id'.")
            }
            else if (row.sample_id.contains(" ")) {
                error("Invalid value in csv file. Value for 'sample_id' can not contain space.")
            }
            [ [ row.sample_id.toString() ], row ]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [ rows, size ]
        }.transpose()
        .map { row, num_lanes ->
            def meta = [:]
            def sample_id = row.sample_id
            meta.sample_id = sample_id
            meta.run_id = run_id
            def sample_type = 'tumor'
            if (row.sample_type) {
                sample_type = row.sample_type
            }
            meta.sample_type = sample_type

            // mapping with fastq
            if (row.fastq_2) {
                def fastq_1 = file(row.fastq_1, checkIfExists: true)
                def fastq_2 = file(row.fastq_2, checkIfExists: true)
                def platform = row.platform ?: "ILLUMINA"
                def (flowcell, lane, machine, run_nr) = flowcellLaneFromFastq(fastq_1)
                def rg_id = "${sample_id}_${flowcell}_${lane}"
                def read_group  = "\"@RG\\tID:${rg_id}\\tSM:${sample_id}\\tPL:${platform}\\tLB:${sample_id}\\tPU:${flowcell}\""
                meta.id = rg_id.toString()                
                meta.read_group = read_group.toString()
                meta.data_type = "fastq"
                return [ meta, [ fastq_1, fastq_2 ] ]
            }
            // bam part 
            if (row.bam) {
                def bam = file(row.bam, checkIfExists: true)
                meta.data_type = "bam"
                meta.id = meta.sample_id
                def bai = file(row.bai, checkIfExists: true)
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
