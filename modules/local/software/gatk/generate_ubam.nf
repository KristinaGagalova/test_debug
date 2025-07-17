process GENERATE_UBAM {
    tag "${meta.id}"
    label 'medium_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1' :
        'broadinstitute/gatk:4.2.6.1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.ubam"), emit: ubam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def platform = task.ext.platform ?: "ILLUMINA"
    def sequencing_center = task.ext.sequencing_center ?: "Unknown"
    def run_date = task.ext.run_date ?: new Date().format("yyyy-MM-dd'T'HH:mm:ss")
    
    """
    mkdir -p tmp

    gatk FastqToSam \\
        --FASTQ ${reads[0]} \\
        --FASTQ2 ${reads[1]} \\
        --OUTPUT ${prefix}.ubam \\
        --READ_GROUP_NAME ${prefix} \\
        --SAMPLE_NAME ${prefix} \\
        --LIBRARY_NAME ${prefix} \\
        --PLATFORM ${platform} \\
        --SEQUENCING_CENTER ${sequencing_center} \\
        --RUN_DATE ${run_date} \\
        --TMP_DIR tmp \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^(GATK|The Genome Analysis Toolkit)' | sed 's/.*v//' | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ubam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}