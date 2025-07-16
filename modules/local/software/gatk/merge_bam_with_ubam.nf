process MERGE_BAM_WITH_UBAM {
    
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
    tuple val(meta), path(reference)
    tuple val(meta), path(seq_dict)

    output:
    path ubam
    path "${bam.simpleName}.final.bam"
    path "${bam.simpleName}.final.bai"

    """
    echo ${ubam}
    mkdir tmp
    gatk MergeBamAlignment \
	    --ALIGNED ${bam} \
        --UNMAPPED ${ubam} \
	    -O ${bam.simpleName}.final.bam \
	    -R ${reference}
    gatk BuildBamIndex \
	    -I ${bam.simpleName}.final.bam \
	    -O ${bam.simpleName}.final.bai
    """
}