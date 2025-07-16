process MARK_DUPLICATES {

    label 'mediumHighMem_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path bam
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${bam.baseName}.marked_duplicates.bam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    gatk MarkDuplicates \
	    I=${bam} \
	    O=${bam.baseName}.marked_duplicates.bam \
	    M=${bam.baseName}.marked_duplicates.txt
    """
}