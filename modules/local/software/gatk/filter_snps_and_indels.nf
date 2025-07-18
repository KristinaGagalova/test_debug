process FILTER_SNPS_AND_INDELS {

    tag "${meta.id}"
    label 'medium_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta1), path(vcf_index)
    tuple val(ref_meta), path(reference)
    tuple val(fai_meta), path(fai_index)
    // tuple val(dict_meta), path(seq_dict)

    output:
    tuple val(meta), path("*.snps.vcf")     , emit: snps_vcf
    tuple val(meta), path("*.snps.vcf.idx") , emit: snps_vcf_index
    tuple val(meta), path("*.indels.vcf")   , emit: indels_vcf
    tuple val(meta), path("*.indels.vcf.idx"), emit: indels_vcf_index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // SelectVariants SNPs args
    def args2 = task.ext.args2 ?: ''  // SelectVariants INDELs args
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def ref_name = params.refname ?: meta.id
    // Use meta to avoid unused variable error
    println "Processing meta: ${meta}"
    """
    mkdir -p ${tmp_dir}

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        --output ${ref_name}.combined_panel.snps.vcf \\
        --tmp-dir ${tmp_dir} \\
        ${args1}

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include INDEL \\
        --output ${ref_name}.combined_panel.indels.vcf \\
        --tmp-dir ${tmp_dir} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def ref_name = params.refname ?: meta.id
    """
    touch ${ref_name}.combined_panel.snps.vcf
    touch ${ref_name}.combined_panel.snps.vcf.idx
    touch ${ref_name}.combined_panel.indels.vcf
    touch ${ref_name}.combined_panel.indels.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}' || echo "4.2.6.1")
    END_VERSIONS
    """
}