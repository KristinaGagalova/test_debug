process GENERATE_UBAM {
    
    label 'medium_task'
    label 'gatk'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"
    
    input:
    tuple val(sampleId), file(reads)

    output:
    tuple val(sampleId), file(reads)
    path "${sampleId}.ubam"

    """
    mkdir tmp
    platform="Illumina"
    seqprovider="AGRF"
    date="2023-01-40T00:00:00-0400"
    gatk FastqToSam \
        -F1 ${reads[0]} \
        -F2 ${reads[1]} \
        -O ${sampleId}".ubam" \
        --READ_GROUP_NAME ${sampleId} \
        --SAMPLE_NAME ${sampleId}\
        --LIBRARY_NAME ${sampleId} \
        --PLATFORM \$platform \
        --SEQUENCING_CENTER \$seqprovider \
        --RUN_DATE \$date \
        --TMP_DIR tmp
    """
}
