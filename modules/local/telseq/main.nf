process TELSEQ {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/telseq:0.0.2--ha7703dc_5':
        'biocontainers/telseq:0.0.2--ha7703dc_5' }"

    input:
    tuple val(meta), path(bam), path(bai)
   
    output:
    path "*.telseq.txt", emit: output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    telseq \\
        ${args} \\
        -o ${prefix}.telseq.txt \\
        ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telseq: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch {prefix}.telseq.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telseq: ${VERSION}
    END_VERSIONS
    """
}