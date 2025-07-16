process QUALITY_FILTER_VARIANTS {
    
    tag "${meta.id}"
    label 'mediumHighMem_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(snp_vcf)
    tuple val(meta), path(snp_vcf_index)
    tuple val(meta), path(indel_vcf)
    tuple val(meta), path(indel_vcf_index)
    tuple val(ref_meta), path(reference)
    tuple val(fai_meta), path(fai_index)
    tuple val(dict_meta), path(seq_dict)

    output:
    tuple val(meta), path("*.filtered.snps.vcf")     , emit: filtered_snps_vcf
    tuple val(meta), path("*.filtered.snps.vcf.idx") , emit: filtered_snps_vcf_index
    tuple val(meta), path("*.filtered.indels.vcf")   , emit: filtered_indels_vcf
    tuple val(meta), path("*.filtered.indels.vcf.idx"), emit: filtered_indels_vcf_index
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // VariantFiltration SNPs args
    def args2 = task.ext.args2 ?: ''  // VariantFiltration INDELs args
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def ref_name = params.refname ?: meta.id
    
    // SNP filter parameters
    def snp_qd = task.ext.snp_qd ?: "2.0"
    def snp_qual = task.ext.snp_qual ?: "30.0"
    def snp_sor = task.ext.snp_sor ?: "3.0"
    def snp_fs = task.ext.snp_fs ?: "60.0"
    def snp_mq = task.ext.snp_mq ?: "40.0"
    def snp_mqranksum = task.ext.snp_mqranksum ?: "-12.5"
    def snp_readposranksum = task.ext.snp_readposranksum ?: "-8.0"
    
    // INDEL filter parameters
    def indel_qd = task.ext.indel_qd ?: "2.0"
    def indel_qual = task.ext.indel_qual ?: "30.0"
    def indel_fs = task.ext.indel_fs ?: "200.0"
    def indel_readposranksum = task.ext.indel_readposranksum ?: "-20.0"
    """
    mkdir -p ${tmp_dir}

    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${snp_vcf} \\
        -filter "QD < ${snp_qd}" --filter-name "QD${snp_qd}" \\
        -filter "QUAL < ${snp_qual}" --filter-name "QUAL${snp_qual}" \\
        -filter "SOR > ${snp_sor}" --filter-name "SOR${snp_sor}" \\
        -filter "FS > ${snp_fs}" --filter-name "FS${snp_fs}" \\
        -filter "MQ < ${snp_mq}" --filter-name "MQ${snp_mq}" \\
        -filter "MQRankSum < ${snp_mqranksum}" --filter-name "MQRankSum${snp_mqranksum}" \\
        -filter "ReadPosRankSum < ${snp_readposranksum}" --filter-name "ReadPosRankSum${snp_readposranksum}" \\
        -O ${ref_name}.combined_panel.filtered.snps.vcf \\
        --tmp-dir ${tmp_dir} \\
        ${args1}
	
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${indel_vcf} \\
        -filter "QD < ${indel_qd}" --filter-name "QD${indel_qd}" \\
        -filter "QUAL < ${indel_qual}" --filter-name "QUAL${indel_qual}" \\
        -filter "FS > ${indel_fs}" --filter-name "FS${indel_fs}" \\
        -filter "ReadPosRankSum < ${indel_readposranksum}" --filter-name "ReadPosRankSum${indel_readposranksum}" \\
        -O ${ref_name}.combined_panel.filtered.indels.vcf \\
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
    touch ${ref_name}.combined_panel.filtered.snps.vcf
    touch ${ref_name}.combined_panel.filtered.snps.vcf.idx
    touch ${ref_name}.combined_panel.filtered.indels.vcf
    touch ${ref_name}.combined_panel.filtered.indels.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}' || echo "4.2.6.1")
    END_VERSIONS
    """
}