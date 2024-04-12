process HMFTOOLS_LINX_SVVISUALISER {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-linx:1.23.6--hdfd78af_0':
        'biocontainers/hmftools-linx:1.23.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(linx_dir)
    path(ensembl_data_dir)
    
    output:
    tuple val(meta), path("plot/"), emit: linx_plot_dir
    tuple val(meta), path("data/"), emit: linx_data_dir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_id = meta.tumor_sample_id ?: meta.id
    def prefix = task.ext.prefix ?: "${sample_id}"    
    def VERSION = '1.23.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """ 
    java -Xmx${task.memory.toGiga() - 1}g -cp /usr/local/share/hmftools-linx-1.23.6-0/linx.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \\
        -sample ${prefix} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -plot_out .//plot \\
        -data_out ./data \\
        -vis_file_dir ${linx_dir} \\
        -circos circos \\
        -threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def sample_id = meta.tumor_sample_id ?: meta.id
    def prefix = task.ext.prefix ?: "${sample_id}"
    def VERSION = '1.23.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir plot
    mkdir data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: ${VERSION}
    END_VERSIONS
    """
}
