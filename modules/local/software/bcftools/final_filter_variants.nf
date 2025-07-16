process FINAL_FILTER_VARIANTS {

    label 'small_task'
    label 'bcftools'

    conda (params.enable_conda ? "bioconda::bcftools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'staphb/bcftools:1.17':
        null }"

    input:
    path snps
    path snps_idx
    path indels
    path indels_idx

    output:
    path "${params.refname}.combined_panel.filtered.biallelic.1random5kb.ld90.snps.vcf"

    """
    bcftools view \
	--max-alleles 2 \
	--exclude-types indels \
	-e "F_MISSING>=${params.missingdata}||FILTER!='PASS'" \
	${snps} -o ${params.refname}.combined_panel.filtered.biallelic.snps.bcf
	
    bcftools +prune \
        -m 0.9 \
        -w 5000bp \
        ${params.refname}.combined_panel.filtered.biallelic.snps.bcf \
        -Ov \
        -o ${params.refname}.combined_panel.filtered.biallelic.1random5kb.ld90.snps.vcf
    """
}