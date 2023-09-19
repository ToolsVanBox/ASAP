process HMFTOOLS_PURPLE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-purple:3.7.1--hdfd78af_0':
        'biocontainers/hmftools-purple:3.7.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(gripss_vcf), path(gripss_vcf_tbi), path(gripss_filtered_vcf), path(gripss_filtered_vcf_tbi), path(amber), path(cobalt)
    path(fasta)
    path(fai)
    path(fasta_dict)
    path(gc_profile)
    path(driver_gene_panel)
    path(known_hotspot_file)
    path(ensembl_data_dir)
    
    output:
    tuple val(meta), path("*driver.catalog.somatic.tsv"), emit: driver_catalog_somatic_tsv
    tuple val(meta), path("*purple.cnv.gene.tsv"), emit: purple_cnv_gene_tsv
    tuple val(meta), path("*purple.cnv.somatic.tsv"), emit: purple_cnv_somatic_tsv
    tuple val(meta), path("*purple.purity.range.tsv"), emit: purple_purity_range_tsv
    tuple val(meta), path("*purple.purity.tsv"), emit: purple_purity_tsv
    tuple val(meta), path("*purple.qc"), emit: purple_qc
    tuple val(meta), path("*purple.segment.tsv"), emit: purple_segment_tsv
    tuple val(meta), path("*purple.somatic.clonality.tsv"), emit: purple_somatic_clonality_tsv
    tuple val(meta), path("*purple.sv.vcf.gz"), emit: purple_sv_vcf
    tuple val(meta), path("*purple.sv.vcf.gz.tbi"), emit: purple_sv_vcf_tbi    
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def gc_profile = gc_profile ? "-gc_profile ${gc_profile}" : ""
    def run_drivers = driver_gene_panel ? "-run_drivers -driver_gene_panel ${driver_gene_panel}" : ""
    def known_hotspot_file = known_hotspot_file ? "-somatic_hotspots ${known_hotspot_file}" : ""
    def VERSION = '3.7.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """ 
    PURPLE \\
        ${args} \\
        -Xmx${task.memory.toGiga() - 1}g \\
        -threads ${task.cpus} \\
        -tumor ${prefix} \\
        -amber ${amber} \\
        -cobalt ${cobalt} \\
        -output_dir ./ \\
        ${gc_profile} \\
        -ensembl_data_dir ${ensembl_data_dir} \\
        -ref_genome ${fasta} \\
        -structural_vcf ${gripss_filtered_vcf} \\
        -sv_recovery_vcf ${gripss_vcf} \\
        ${run_drivers} \\
        ${known_hotspot_file} \\
        -circos circos
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '3.7.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.driver.catalog.somatic.tsv
    touch ${prefix}.purple.cnv.gene.tsv
    touch ${prefix}.purple.cnv.somatic.tsv
    touch ${prefix}.purple.purity.range.tsv
    touch ${prefix}.purple.qc
    touch ${prefix}.purple.segment.tsv
    touch ${prefix}.purple.somatic.clonality.tsv
    touch ${prefix}.purple.sv.vcf.gz
    touch ${prefix}.purple.sv.vcf.gz.tbi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: ${VERSION}
    END_VERSIONS
    """
}
