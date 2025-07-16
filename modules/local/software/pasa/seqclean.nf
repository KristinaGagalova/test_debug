//reference code: https://github.com/nf-core/genomeannotator

process PASA_SEQCLEAN {

    tag "${meta}"
    label 'pasa'
    label 'big_task'
    
    container 'pasapipeline/pasapipeline:2.5.2'
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("*.clean"), path("*.cln"), emit: fasta
    path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    // Seqclean won't work inside container if the username is not exported to USER
    """
    export USER=${workflow.userName}
    seqclean $fasta -c 6 #${task.cpus} hardocoded CPUs - launches an error with more CPUs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
		pasa: 2.5.2
    END_VERSIONS
    """
}
