process DICT_REF {

    label 'small_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(meta), path(reference)  

    output:
    tuple val(meta), path("*.dict"), emit:index

    """
    gatk CreateSequenceDictionary -R ${reference}
    """
}