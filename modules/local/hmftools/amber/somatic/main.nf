process HMFTOOLS_AMBER_SOMATIC {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-amber:3.9--hdfd78af_1':
        'biocontainers/hmftools-amber:3.9--hdfd78af_1' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path(fasta)
    path(bafsnps)
    
    output:
    tuple val(meta), path("*amber.baf.pcf") , emit: amber_baf_pcf
    tuple val(meta), path("*amber.baf.tsv.gz") , emit: amber_baf_tsv_gz
    tuple val(meta), path("*amber.contamination.tsv") , emit: amber_contamination_tsv
    tuple val(meta), path("*amber.contamination.vcf.gz") , emit: amber_contamination_vcf_gz
    tuple val(meta), path("*amber.contamination.vcf.gz.tbi") , emit: amber_contamination_vcf_gz_tbi
    tuple val(meta), path("*amber.qc") , emit: amber_qc
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"    
    def args = task.ext.args ?: ''
    def tumor = meta.tumor_sample_id
    def reference = meta.normal_sample_id
    bafsnps = bafsnps ? "-loci ${bafsnps}" : ""
    def VERSION = '3.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    AMBER \\
        ${args} \\
        -Xmx${task.memory.toGiga() - 1}g \\
		-threads ${task.cpus} \\
		-tumor ${tumor} \\
        -tumor_bam ${tumor_bam} \\
        -reference ${reference} \\
        -reference_bam ${normal_bam} \\
		${bafsnps} \\
		-ref_genome ${fasta} \\
		-output_dir ./
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor = meta.tumor_sample_id    
    def args = task.ext.args ?: ''    
    def VERSION = '3.9' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${tumor}.amber.baf.pcf
    touch ${tumor}.amber.baf.tsv.gz
    touch ${tumor}.amber.qc
    touch ${tumor}.amber.contamination.tsv
    touch ${tumor}.amber.contamination.vcf.gz
    touch ${tumor}.amber.contamination.vcf.gz.tbi
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amber: ${VERSION}
    END_VERSIONS
    """
}
