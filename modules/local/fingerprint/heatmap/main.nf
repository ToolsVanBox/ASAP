process FINGERPRINT_HEATMAP {
    tag "$meta.id"
    label 'process_single'
    // container = 'docker.io/vanboxtelbioinformatics/asap_r:1.0'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/asap_r:1.1':
        params.artifact_registry_path + '/asap_r@sha256:17e59124425ae3aa17fe096418f1d61f225e9c3cae163430fbeafbc8d716d260' }"


    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*_fingerprintheatmap.pdf")     , emit: fingerprintheatmap
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input_vcfs = vcf.join(',')
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R --slave --file=/ASAP_R/fingerprintheatmap.R --args ${input_vcfs} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) |  grep -oP "R version .+ --" | cut -f 3 -d' ' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_fingerprintheatmap.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	R: \$(echo \$(R --version 2>&1) |  grep -oP "R version .+ --" | cut -f 3 -d' ' )
    END_VERSIONS
    """

}