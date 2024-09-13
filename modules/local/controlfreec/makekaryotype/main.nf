process CONTROLFREEC_MAKEKARYOTYPE {
    tag "$meta.id"
    label 'process_low'
    container = 'docker://vanboxtelbioinformatics/asap_r:1.0'

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*karyotype.pdf")     , emit: karyotype
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ploidy = task.ext.ploidy ?: 2
    def maxLevelToPlot = task.ext.maxleveltoplot ?: 4
    def binsize = task.ext.binsize ?: 50000
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R --slave --file=/ASAP_R/plotkaryotype.R --args ${ploidy} ${maxLevelToPlot} ${binsize} ${ratio} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) ) |  grep -oP "R version .+ --" | cut -f 3 -d' '
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_karyotype.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) ) |  grep -oP "R version .+ --" | cut -f 3 -d' '
    END_VERSIONS
    """
}
