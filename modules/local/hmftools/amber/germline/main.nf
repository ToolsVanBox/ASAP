process HMFTOOLS_AMBER_GERMLINE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:3.9--hdfd78af_1':
        'biocontainers/hmftools-amber:3.9--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(bafsnps)
    
    output:
    tuple val(meta), path("*amber.homozygousregion.tsv") , emit: amber_homozygousregion_tsv
    tuple val(meta), path("*amber.baf.tsv.gz") , emit: amber_baf_tsv_gz
    tuple val(meta), path("*amber.qc") , emit: amber_qc
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    def reference = meta.id
    def bafsnps = bafsnps ? "-loci ${bafsnps}" : ""
    def VERSION = '3.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    AMBER \\
        ${args} \\
        -Xmx${task.memory.toGiga() - 1}g \\
		-threads ${task.cpus} \\
		-reference ${reference} \\
        -reference_bam ${bam} \\
        ${bafsnps} \\
		-ref_genome ${fasta} \\
		-output_dir ./
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"    
    
    def VERSION = '3.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.amber.baf.pcf
    touch ${prefix}.amber.baf.tsv.gz
    touch ${prefix}.amber.qc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: ${VERSION}
    END_VERSIONS
    """
}
