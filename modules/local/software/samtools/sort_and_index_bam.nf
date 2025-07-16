process SORT_AND_INDEX_BAM {
    
    label 'small_task'
    label 'samtools'

    conda (params.enable_conda ? "bioconda::samtools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/samtools:1.17--hd87286a_2':
        null }"

    input:
    path sam
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${sam.baseName}_v_${params.refname}.bam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    samtools sort ${sam} \
		-O BAM | \
		tee ${sam.baseName}_v_${params.refname}.bam | \
		samtools index - ${sam.baseName}_v_${params.refname}.bam.bai
    """
}