process CONTROLFREEC_MAKEBAFPLOT {
    tag "$meta.id"
    label 'process_low'
    // container = 'docker.io/vanboxtelbioinformatics/asap_r:1.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/asap_r:1.1':
        params.artifact_registry_path + '/asap_r:1.1' }"
    input:
    tuple val(meta), path(baf)

    output:
    tuple val(meta), path("*_BAF.pdf")          , emit: bafplot
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def binsize = task.ext.binsize ?: 1000
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R --slave --file=/ASAP_R/plotbaf.R --args ${baf} ${binsize} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) |  grep -oP "R version .+ --" | cut -f 3 -d' ' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_BAF.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	R: \$(echo \$(R --version 2>&1) |  grep -oP "R version .+ --" | cut -f 3 -d' ' )
    END_VERSIONS
    """
}
