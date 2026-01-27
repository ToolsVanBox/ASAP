process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    // conda "bioconda::gatk4=4.4.0.0 conda-forge::openjdk=8.0.312"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0':
        params.artifact_registry_path + '/gatk4@sha256:573f43fe59145d65106e5b2647d5c0f1662820716234a04104dfe1c91ebfa4be' }"
//        'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' }"
        
    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fasta_fai
    path  dict

    output:
    tuple val(meta), path("${prefix}.bam"),     emit: output
    tuple val(meta), path("${prefix}.bai"), emit: bam_index, optional:true
    tuple val(meta), path("*.metrics"),     emit: metrics, optional: true
    path "versions.yml"               ,     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bam.collect{"--input $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" MarkDuplicatesSpark \\
        $input_list \\
        --output ${prefix}.bam \\
        --reference $fasta \\
        --spark-master local[${task.cpus}] \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        openjdk: \$(echo \$(java -version 2>&1) | grep version | sed 's/\"//g' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
