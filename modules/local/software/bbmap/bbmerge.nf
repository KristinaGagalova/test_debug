process BBMERGE {

    label 'small_task'
    label 'bbmap'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::bbmap=38.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bbmap:38.19--h470a237_0':
        null }"

    input: 
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("unmerged_${reads1}"), path("unmerged_${reads2}"), path("merged_${reads1}"), emit: reads
    tuple val(meta), path("${meta}_ihist.txt"), emit: insert_histogram

    """
    bbmerge.sh \\
        --in1=${reads1} \\
        --in2=${reads2} \\
        out=merged_${reads1} \\
        outu1=unmerged_${reads1} \\
        outu2=unmerged_${reads2} \\
        ihist=${meta}_ihist.txt
    """

}