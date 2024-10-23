process SMURF {
  tag "$meta.id"
  label 'process_single'
  
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/smurf:3.0.2':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/smurf@sha256:3ea5c21b0aaf696dbdaadc851ada4fc82744ece9b03f139bca73a7299e215492' }"
  // container = 'docker.io/vanboxtelbioinformatics/smurf:3.0.2'

  input:
    tuple val(meta), path(vcf), path(tbi), path(bams), path(bai)    
    path( config )

  output:
    tuple val(meta), path("*SMuRF.vcf"), emit: smurf_vcf
    tuple val(meta), path("*SMuRF.filtered.vcf"), emit: smurf_filtered_vcf
    tuple val(meta), path("*.pdf"), emit: smurf_pdf
    path "versions.yml", emit: versions


  script:
    b = bams ? ' -b ' + bams.join(' -b ') : ''
    n = meta.bulk_names ? ' -n ' + meta.bulk_names.join(' -n ') : ''
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def config = ${params.genomes[params.genome].smurf_config}
    
    """
    host=\$(hostname)
    echo \${host}

    export projectDir=${projectDir}

    #python /smurf/SMuRF.py \
    python /hpc/pmc_vanboxtel/tools/ToolsVanBox_DEV/SMuRF/SMuRF.py \
    -i ${vcf} \
    ${b} \
    ${n} \
    -t ${task.cpus} \
    -c ${config}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SMuRF: \$(echo \$(python /smurf/SMuRF.py --version 2>&1) )
    END_VERSIONS
    """
}
