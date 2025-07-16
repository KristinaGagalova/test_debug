process VARIANT_CALLING {

    label 'very_large_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path ubam
    path bam
    path bamindex
    tuple val(meta), path(reference)
    tuple val(meta), path(fai_index)
    tuple val(meta), path(seq_dict)

    output:
    path "${bam.simpleName}.g.vcf" 

    """
    mkdir tmp
    gatk HaplotypeCaller \
        -I ${bam} \
        -O ${bam.simpleName}.g.vcf \
        -R ${reference} \
        -ERC GVCF \
        --minimum-mapping-quality 20 \
        --min-base-quality-score 20 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        --tmp-dir tmp \
        --native-pair-hmm-threads ${task.cpus * 2}
    """ 
}