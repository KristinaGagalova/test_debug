process SAMTOOLS_INDEX_GENOME {

    label 'samtools'
    label 'small_task'

    conda (params.enable_conda ? "bioconda::samtools=1.21" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0':
        null }"
        
    input:
    tuple val(meta), path(reference) 

    output:
    tuple val(meta), path("${reference}.fai"), emit:index

    """
    samtools faidx ${reference}
    """
}