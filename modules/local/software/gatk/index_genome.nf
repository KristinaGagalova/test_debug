process DICT_REF {

    tag "${meta.id}"
    label 'small_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(reference)  

    output:
    tuple val(meta), path("${meta.id}.dict"), emit: index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk CreateSequenceDictionary \\
        -R ${reference} \\
        -O ${prefix}.dict \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${reference.baseName}"
    """
    touch ${prefix}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """
}
