process BBSPLIT {

    label 'large_task'
    label 'bbmap'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::bbmap=38.19" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bbmap:38.19--h470a237_0':
        null }"

    input: 
    tuple val(meta), path(cleaned_reads1), path(cleaned_reads2)
    path(reference)
    val(bin)

    output:
    tuple val(meta), path("${meta}_${bin}_1.fq"), path("${meta}_${bin}_2.fq"), emit: reads_in_bin

    script:
    def all_references = reference.collect { it.toString() }.join(",") 

    """
    echo ${all_references}

    bbsplit.sh \\
        ref=${all_references} \\
        nzo=f \\
        in1=${cleaned_reads1}  \\
        in2=${cleaned_reads2} \\
        outu1=${meta}_1.unknown.fastq \\
        outu2=${meta}_2.unknown.fastq \\
        basename=${meta}_%_#.fq 
    """

}