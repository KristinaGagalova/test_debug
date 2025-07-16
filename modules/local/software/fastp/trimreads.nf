process TRIM_READS {

    label 'fastp'
    label 'small_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::fastp=0.23.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0':
        null }"

    input: 
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("trimmed_${reads1}"), path("trimmed_${reads2}")

    """
    fastp \\
        --in1 ${reads1} \\
        --in2 ${reads2} \\
        --out1 trimmed_${reads1} \\
        --out2 trimmed_${reads2}
    """

}