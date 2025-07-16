process MARK_ILLUMINA_ADAPTERS {

    label 'small_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${sampleId}.marked.fastq"
    tuple val(sampleId), file(reads)
    path ubam

    """
    mkdir tmp
    gatk MarkIlluminaAdapters \
        -I ${ubam} \
        -O ${sampleId}.marked.ubam \
        -M ${sampleId}.markilluminaadapters_metrics.txt \
        -TMP_DIR tmp
    gatk SamToFastq \
		-I ${sampleId}.marked.ubam \
		-FASTQ ${sampleId}.marked.fastq \
		-CLIPPING_ATTRIBUTE XT \
		-CLIPPING_ACTION 2 \
		-INTERLEAVE true \
		-NON_PF true 
    """
}