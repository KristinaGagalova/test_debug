process COMBINE_AND_GENOTYPE_VCF {

    label 'mediumHighMem_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path gvcfs
    tuple val(meta), path(reference)
    tuple val(meta), path(fai_index)
    tuple val(meta), path(seq_dict)

    output:
    path "${params.refname}.combined_panel.vcf"
    path "${params.refname}.combined_panel.vcf.idx"

    """
    ls *.g.vcf > list.gvcfs.list
    gatk CombineGVCFs \
	    -R ${reference} \
	    --variant list.gvcfs.list \
	    -O ${params.refname}.combined_panel.g.vcf
    gatk GenotypeGVCFs \
	    -R ${reference} \
	    --variant ${params.refname}.combined_panel.g.vcf \
	    -O ${params.refname}.combined_panel.vcf
    """
}
