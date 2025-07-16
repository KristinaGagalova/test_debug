process FILTER_SNPS_AND_INDELS {

    label 'medium_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path vcf
    path vcf_index
    tuple val(meta), path(reference)
    tuple val(meta), path(fai_index)
    tuple val(meta), path(seq_dict)

    output:
    path "${params.refname}.combined_panel.snps.vcf"
    path "${params.refname}.combined_panel.snps.vcf.idx"
    path "${params.refname}.combined_panel.indels.vcf"
    path "${params.refname}.combined_panel.indels.vcf.idx"

    """
    gatk SelectVariants \
        -R ${reference} \
        -V ${vcf} \
        --select-type-to-include SNP \
        --output ${params.refname}.combined_panel.snps.vcf
    gatk SelectVariants \
        -R ${reference} \
        -V ${vcf} \
        --select-type-to-include INDEL \
        --output ${params.refname}.combined_panel.indels.vcf
    """
}
