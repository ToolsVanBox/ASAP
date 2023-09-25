process HMFTOOLS_LINX_TUMORONLY {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-linx:1.23.6--hdfd78af_0':
        'biocontainers/hmftools-linx:1.23.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(purple_sv_vcf), path(purple_sv_vcf_tbi), path(purple_dir)
    path(ensembl_data_dir)
    path(driver_gene_panel)
    path(known_fusion_file)
    
    output:
    tuple val(meta), path("*linx.breakend.tsv"), emit: linx_breakend_tsv
    tuple val(meta), path("*linx.clusters.tsv"), emit: linx_clusters_tsv
    tuple val(meta), path("*linx.driver.catalog.tsv"), emit: linx_driver_catalog_tsv
    tuple val(meta), path("*linx.drivers.tsv"), emit: linx_drivers_tsv
    tuple val(meta), path("*linx.fusion.tsv"), emit: linx_fusion_tsv
    tuple val(meta), path("*linx.links.tsv"), emit: linx_links_tsv
    tuple val(meta), path("*linx.svs.tsv"), emit: linx_svs_tsv
    tuple val(meta), path("*linx.vis_copy_number.tsv"), emit: linx_vis_copy_number_tsv
    tuple val(meta), path("*linx.vis_fusion.tsv"), emit: linx_vis_fusion_tsv
    tuple val(meta), path("*linx.vis_gene_exon.tsv"), emit: linx_vis_gene_exon_tsv
    tuple val(meta), path("*linx.vis_protein_domain.tsv"), emit: linx_vis_protein_domain_tsv
    tuple val(meta), path("*linx.vis_segments.tsv"), emit: linx_vis_segments_tsv
    tuple val(meta), path("*linx.vis_sv_data.tsv"), emit: linx_vis_sv_data_tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def sample = meta.id
    def driver_gene_panel = driver_gene_panel ? "-driver_gene_panel ${driver_gene_panel}" : ""
    def known_fusion_file = known_fusion_file ? "-known_fusion_file ${known_fusion_file}" : ""
    def VERSION = '1.23.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """ 
    linx \\
        ${args} \\
        -Xmx${task.memory.toGiga() - 1}g \\
        -sample ${sample} \\
        -sv_vcf ${purple_sv_vcf} \\
        -purple_dir ${purple_dir} \\
        -output_dir ./ \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        ${known_fusion_file} \\
        ${driver_gene_panel}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def VERSION = '1.23.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.linx.breakend.tsv
    touch ${prefix}.linx.clusters.tsv
    touch ${prefix}.linx.driver.catalog.tsv
    touch ${prefix}.linx.drivers.tsv
    touch ${prefix}.linx.fusion.tsv
    touch ${prefix}.linx.links.tsv
    touch ${prefix}.linx.svs.tsv
    touch ${prefix}.linx.vis_copy_number.tsv
    touch ${prefix}.linx.vis_fusion.tsv
    touch ${prefix}.linx.vis_gene_exon.tsv
    touch ${prefix}.linx.vis_protein_domain.tsv
    touch ${prefix}.linx.vis_segments.tsv
    touch ${prefix}.linx.vis_sv_data.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        linx: ${VERSION}
    END_VERSIONS
    """
}
