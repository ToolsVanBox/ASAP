
nextflow.enable.dsl=2

include { FASTQ_QC_PRE_MAPPING } from '../subworkflows/local/fastq_qc_pre_mapping/main.nf'
include { FASTQ_ALIGN } from '../subworkflows/local/fastq_align/main.nf'

input_sample = parseSampleSheet( params.input )
input_sample_type = input_sample.branch{
        bam:   it[0].data_type == "bam"
        fastq: it[0].data_type == "fastq"
    }

workflow WAP {
    // FASTQ_QC_PRE_MAPPING( input_sample_type.fastq )    

    FASTQ_ALIGN( input_sample_type.fastq )

}

def parseSampleSheet(csv_file) {
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
                def bai = file(row.bai, checkIfExists: true)
                return [ meta, [bam, bai] ]
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
