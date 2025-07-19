process FILTER_SNPS_AND_INDELS {

    tag "${meta.id}"
    label 'medium_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path(reference)
    path(fai_index)
    path(seq_dict)

    output:
    tuple val(meta),
	path("${meta.id}.snps.vcf")   , path("${meta.id}.snps.vcf.idx"),
	path("${meta.id}.indels.vcf") , path("${meta.id}.indels.vcf.idx") , emit: vcfs
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // SelectVariants SNPs args
    def args2 = task.ext.args2 ?: ''  // SelectVariants INDELs args
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def prefix = meta.id
    """
    mkdir -p ${tmp_dir}

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        --output ${prefix}.snps.vcf \\
        --tmp-dir ${tmp_dir} \\
        ${args1}

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        --output ${prefix}.indels.vcf \\
        --tmp-dir ${tmp_dir} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}.snps.vcf
    touch ${prefix}.snps.vcf.idx
    touch ${prefix}.indels.vcf
    touch ${prefix}.indels.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """
}
