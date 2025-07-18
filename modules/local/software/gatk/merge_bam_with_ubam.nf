process MERGE_BAM_WITH_UBAM {
    
    tag "${meta.id}"
    label 'mediumHighMem_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta3), path(ubam)
    tuple val(ref_meta), path(reference)

    output:
    tuple val(meta), path("*.final.bam"), emit: bam
    tuple val(meta), path("*.final.bai"), emit: bai
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // MergeBamAlignment args
    def args2 = task.ext.args2 ?: ''  // BuildBamIndex args
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def sort_order = task.ext.sort_order ?: "coordinate"
    def clip_adapters = task.ext.clip_adapters ? "true" : "false"
    """
    echo "Processing: ${ubam}"
    mkdir -p ${tmp_dir}

    gatk MergeBamAlignment \\
        --ALIGNED ${bam} \\
        --UNMAPPED ${ubam} \\
        -O ${prefix}.final.bam \\
        -R ${reference} \\
        --TMP_DIR ${tmp_dir} \\
        --SO ${sort_order} \\
        --CREATE_INDEX false \\
        --CLIP_ADAPTERS ${clip_adapters} \\
        ${args1}

    gatk BuildBamIndex \\
        -I ${prefix}.final.bam \\
        -O ${prefix}.final.bai \\
        --TMP_DIR ${tmp_dir} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.final.bam
    touch ${prefix}.final.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}' || echo "4.2.6.1")
    END_VERSIONS
    """
}