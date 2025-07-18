process MARK_DUPLICATES {

    tag "${meta.id}"
    label 'mediumHighMem_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.marked_duplicates.bam")     , emit: bam
    tuple val(meta), path("*.marked_duplicates.txt")     , emit: metrics
    tuple val(meta), path("*.marked_duplicates.bai"), optional: true, emit: bai
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def assume_sort_order = task.ext.assume_sort_order ?: "coordinate"
    def create_index = task.ext.create_index ? "true" : "false"
    """
    mkdir -p ${tmp_dir}

    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${prefix}.marked_duplicates.bam \\
        -M ${prefix}.marked_duplicates.txt \\
        --ASSUME_SORT_ORDER ${assume_sort_order} \\
        --TMP_DIR ${tmp_dir} \\
        --CREATE_INDEX ${create_index} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def create_index = task.ext.create_index ? "true" : "false"
    """
    touch ${prefix}.marked_duplicates.bam
    touch ${prefix}.marked_duplicates.txt
    if [ "${create_index}" == "true" ]; then
        touch ${prefix}.marked_duplicates.bai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6} )
    END_VERSIONS
    """
}
