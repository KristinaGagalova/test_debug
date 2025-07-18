process MARK_ILLUMINA_ADAPTERS {

    tag "${meta.id}"
    label 'small_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(ubam)

    output:
    tuple val(meta), path("${meta.id}.marked.fastq")           , emit: marked_fastq
    tuple val(meta), path(reads)                               , emit: reads
    tuple val(meta), path("${meta.id}.marked.ubam")            , emit: marked_ubam
    tuple val(meta), path("${meta.id}.markilluminaadapters_metrics.txt"), emit: metrics
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def clipping_attribute = task.ext.clipping_attribute ?: "XT"
    def clipping_action = task.ext.clipping_action ?: "2"
    def interleave = task.ext.interleave ? "true" : "false"
    def non_pf = task.ext.non_pf ? "true" : "false"
    """
    mkdir -p ${tmp_dir}
    
    gatk MarkIlluminaAdapters \\
        -I ${ubam} \\
        -O ${prefix}.marked.ubam \\
        -M ${prefix}.markilluminaadapters_metrics.txt \\
        -TMP_DIR ${tmp_dir} \\
        ${args1}
    
    gatk SamToFastq \\
        -I ${prefix}.marked.ubam \\
        -FASTQ ${prefix}.marked.fastq \\
        -CLIPPING_ATTRIBUTE ${clipping_attribute} \\
        -CLIPPING_ACTION ${clipping_action} \\
        -INTERLEAVE ${interleave} \\
        -NON_PF ${non_pf} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.marked.fastq
    touch ${prefix}.markilluminaadapters_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}' || echo "4.2.6.1")
    END_VERSIONS
    """
}