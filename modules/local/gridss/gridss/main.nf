process GRIDSS {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta), path(gridss_bams)
    tuple val(meta_bai), path(gridss_bai)
    tuple val(meta_labels), val(gridss_labels)
    tuple val(meta2), path(bwa_index)
    path(fai)
    path(fasta_dict)
    path(blacklist)
    path(gridss_properties)

    output:
    tuple val(meta), path("*gridss.driver.vcf.gz") , emit: gridss_driver_vcf
    tuple val(meta), path("*gridss.driver.vcf.gz.tbi") , emit: gridss_driver_vcf_tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def bams = gridss_bams.collect{"$it"}.join(' ')
    def labels = gridss_labels.collect{"$it"}.join(',')
    blacklist = blacklist ? "--blacklist ${blacklist}" : ""
    gridss_properties = gridss_properties ? "--configuration ${gridss_properties}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    FASTA=`find -L ./ -name "*.fasta"`

    gridss \\
        ${args} \\
        --output ${prefix}.gridss.driver.vcf.gz \\
        --reference \$FASTA \\
        --threads ${task.cpus} \\
        --assembly ${prefix}.assembly.bam \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        --labels ${labels} \\
        ${blacklist} \\
        ${gridss_properties} \\
        ${bams}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    meta = [ id: ${prefix} ]
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def steps = args.contains("-s ") ? args.split('-s ')[-1].split(" ")[0] :
                args.contains("--steps ") ? args.split('--steps ')[-1].split(" ")[0] :
                "all"
    def vcf = steps.contains("call") || steps.contains("all") ? "touch ${prefix}.vcf.gz" : ""
    def assembly_bam = steps.contains("assembly") || steps.contains("all") ? "touch ${prefix}.assembly.bam" : ""
    """
    ${vcf}
    ${assembly_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
