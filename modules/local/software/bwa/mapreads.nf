process BWA_MAPREADS {

    tag "${meta.id}"  // Fixed to use meta.id
    label 'small_task'
    label 'bwa'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1':
        null }"

    input:
    tuple val(meta), path(reads)  // Can be single interleaved file or paired files
    path(index)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.ext.threads ?: task.cpus
    
    // Handle single interleaved file vs paired files
    def read_input = reads instanceof List ? "${reads[0]} ${reads[1]}" : "${reads}"
    
    """
    bwa mem \\
        -t ${threads} \\
        ${args} \\
        ${reference} \\
        ${read_input} > ${prefix}.sam

    samtools sort \\
        -@ ${threads} \\
        -o "${prefix}.bam" \\
        ${prefix}.sam

    rm ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep "Version" | awk '{print \$2}')
        samtools: \$(samtools version | head -n 1 | cut -d" " -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep "Version" | awk '{print \$2}')
        samtools: \$(samtools version | head -n 1 | cut -d" " -f2)
    END_VERSIONS
    """
}
