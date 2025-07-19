process SAMTOOLS_INDEX_GENOME {

    tag "${meta.id}"
    label 'samtools'
    label 'small_task'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0':
        null }"
        
    input:
    tuple val(meta), path(reference) 

    output:
    tuple val(meta), path("${reference.name}.fai"), emit: index
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools faidx ${args} ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${reference.name}.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | awk '{print \$2}')
    END_VERSIONS
    """
}
