process SEQKIT_DEDUP_FASTA {

    tag "${meta}"
    label "seqkit"
    label 'small_task'

    conda (params.enable_conda ? "bioconda::seqkit=2.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta_dedup), emit: fasta
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    fasta_dedup = fasta.getBaseName() + ".dedup.fasta"
    """
    seqkit rmdup -s < $fasta > $fasta_dedup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        removeduplicated: \$(echo \$(removeduplicated -h | head -n1 | cut -f5 -d " " ))
    END_VERSIONS
    """
}
