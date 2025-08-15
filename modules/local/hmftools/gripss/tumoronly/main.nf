process HMFTOOLS_GRIPSS_TUMORONLY {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-gripss:2.3.2--hdfd78af_0':
        'biocontainers/hmftools-gripss:2.3.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(gridss_unfiltered_vcf)
    tuple val(meta2), path(bwa_index)
    path(fai)
    path(known_hotspot_file)
    path(pon_sgl_file)
	path(pon_sv_file)
    path(repeat_mask_file)

    output:
    tuple val(meta), path("*gripss.vcf.gz") , emit: gripss_vcf
    tuple val(meta), path("*gripss.vcf.gz.tbi") , emit: gripss_vcf_tbi
    tuple val(meta), path("*gripss.filtered.vcf.gz") , emit: gripss_filtered_vcf
    tuple val(meta), path("*gripss.filtered.vcf.gz.tbi") , emit: gripss_filtered_vcf_tbi
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample = meta.id
    def args = task.ext.args ?: ''
    known_hotspot_file = known_hotspot_file ? "-known_hotspot_file ${known_hotspot_file}" : ""
    pon_sgl_file = pon_sgl_file ? "-pon_sgl_file ${pon_sgl_file}" : ""
    pon_sv_file = pon_sv_file ? "-pon_sv_file ${pon_sv_file}" : ""
    repeat_mask_file = repeat_mask_file ? "-repeat_mask_file ${repeat_mask_file}" : ""
    def VERSION = '2.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """    
    FASTA=`find -L ./ -name "*.fasta"`

    gripss \\
        -Xmx${task.memory.toGiga() - 1}g \\
        ${args} \\
        -ref_genome \$FASTA \\
		${known_hotspot_file} \\
		${pon_sgl_file} \\
		${pon_sv_file} \\
        ${repeat_mask_file} \\
		-vcf ${gridss_unfiltered_vcf} \\
		-sample ${sample} \\
        -output_dir ./
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def sample = meta.id
    def VERSION = '2.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${sample}.gripps.vcf.gz
    touch ${sample}.gripps.vcf.gz.tbi
    touch ${sample}.gripps.filtered.vcf.gz
    touch ${sample}.gripps.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: ${VERSION}
    END_VERSIONS
    """
}
