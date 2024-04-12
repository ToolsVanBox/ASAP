process HMFTOOLS_COBALT_TUMORONLY {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-cobalt:1.13--hdfd78af_1':
        'biocontainers/hmftools-cobalt:1.13--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(gc_profile)
    path(tumor_only_diploid_bed)
    
    output:
    tuple val(meta), path("*cobalt.gc.median.tsv") , emit: cobalt_gc_median_tsv
    tuple val(meta), path("*cobalt.ratio.pcf") , emit: cobalt_ratio_pcf
    tuple val(meta), path("*cobalt.ratio.tsv.gz") , emit: cobalt_ratio_tsv_gz
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def tumor = meta.id
    def gc_profile = gc_profile ? "-gc_profile ${gc_profile}" : ""
    def tumor_only_diploid_bed = tumor_only_diploid_bed ? "-tumor_only_diploid_bed ${tumor_only_diploid_bed}" : ""

    def VERSION = '1.13' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """ 
    COBALT \\
        ${args} \\
        -Xmx${task.memory.toGiga() - 1}g \\
        -tumor ${tumor} \\
        -tumor_bam ${bam} \\
        -output_dir ./ \\
        -threads ${task.cpus} \\
        ${gc_profile} \\
        ${tumor_only_diploid_bed}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def VERSION = '1.13' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.cobalt.qc.median.tsv
    touch ${prefix}.cobalt.ratio.pcf
    touch ${prefix}.cobalt.ratio.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobalt: ${VERSION}
    END_VERSIONS
    """
}
