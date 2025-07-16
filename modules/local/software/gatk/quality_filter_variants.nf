process QUALITY_FILTER_VARIANTS {
    
    label 'mediumHighMem_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path snp_vcf
    path snp_vcf_index
    path indel_vcf
    path indel_vcf_index
    tuple val(meta), path(reference)
    tuple val(meta), path(fai_index)
    tuple val(meta), path(seq_dict)

    output:
    path "${params.refname}.combined_panel.filtered.snps.vcf"
    path "${params.refname}.combined_panel.filtered.snps.vcf.idx"
    path "${params.refname}.combined_panel.filtered.indels.vcf"
    path "${params.refname}.combined_panel.filtered.indels.vcf.idx"

    """
    gatk VariantFiltration \
        -R ${reference} \
        -V ${snp_vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${params.refname}.combined_panel.filtered.snps.vcf
	
    gatk VariantFiltration \
        -R ${reference} \
        -V ${indel_vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O ${params.refname}.combined_panel.filtered.indels.vcf
    """
}