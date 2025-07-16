process BWA_INDEX {

    tag "${meta.id}"
    label 'small_task'
    label 'bwa'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1':
        null }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}.*")  , emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bwa index ${args} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         bwa: \$(bwa 2>&1 | grep "Version" | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt
    touch ${fasta}.pac
    touch ${fasta}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         bwa: \$(bwa 2>&1 | grep "Version" | awk '{print \$2}' || echo "0.7.17")
    END_VERSIONS
    """
}