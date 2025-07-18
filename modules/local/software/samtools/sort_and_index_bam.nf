process SORT_AND_INDEX_BAM {
    
    tag "${meta.id}"
    label 'small_task'
    label 'samtools'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/samtools:1.17--hd87286a_2':
        null }"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}*.bam")     , emit: bam
    tuple val(meta), path("${meta.id}*.bam.bai") , emit: bai
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def index_args = task.ext.index_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_suffix = task.ext.ref_suffix ?: (params.refname ? "_v_${params.refname}" : "")
    def output_name = "${prefix}${ref_suffix}.bam"
    def threads = task.ext.threads ?: task.cpus
    """
    samtools sort \\
        -@ ${threads} \\
        -O BAM \\
        ${sam} | \\
        ${args} \\
    tee ${output_name} | \\
    samtools index \\
        ${index_args} \\
        - ${output_name}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -n 1 | cut -d" " -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_suffix = task.ext.ref_suffix ?: (params.refname ? "_v_${params.refname}" : "")
    def output_name = "${prefix}${ref_suffix}.bam"
    """
    touch ${output_name}
    touch ${output_name}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -n 1 | cut -d" " -f2 || echo "1.17")
    END_VERSIONS
    """
}
