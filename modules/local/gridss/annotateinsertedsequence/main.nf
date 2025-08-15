process GRIDSS_ANNOTATEINSERTEDSEQUENCE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta), path(gridss_driver_vcf)
    tuple val(meta2), path(viralreferenc_index)
    path(repeatmaskerbed)

    output:
    tuple val(meta), path("*gridss.unfiltered.vcf.gz")          , emit: gridss_unfiltered_vcf
    tuple val(meta), path("*gridss.unfiltered.vcf.gz.tbi")          , emit: gridss_unfiltered_vcf_tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    repeatmaskerbed = repeatmaskerbed ? "REPEAT_MASKER_BED=${repeatmaskerbed}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """ 
    FASTA=`find -L ./ -name "*.fa"`

    AnnotateInsertedSequence \\
        ${args} \\
        REFERENCE_SEQUENCE=\$FASTA \\
		INPUT=${gridss_driver_vcf} \\
		OUTPUT=${prefix}.gridss.unfiltered.vcf.gz \\
        WORKER_THREADS=${task.cpus} \\
        ${repeatmaskerbed}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_annotateinsertedsequence: ${VERSION}
    END_VERSIONS
    """    

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.gridss.unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_annotateinsertedsequence: ${VERSION}
    END_VERSIONS
    """
}
