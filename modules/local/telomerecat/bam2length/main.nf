process TELOMERECAT_BAM2LENGTH {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/wtsicgp/telomerecat:4.0.2"

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:

    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def VERSION = '4.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    telomerecat bam2length \\
        ${args} \\
        -p ${task.cpus} \\
        --temp_dir ./ \\
        --output ${prefix} \\
        ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerecat: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch blaat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerecat: ${VERSION}
    END_VERSIONS
    """
}
