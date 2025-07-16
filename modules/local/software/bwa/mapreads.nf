process BWA_MAPREADS {

    tag "${meta}"
    label 'small_task'
    label 'bwa'

    conda (params.enable_conda ? "bioconda::bwakit=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1':
        null }"

    input:
    tuple val(meta), path(index)
    tuple val(meta), path(read1), path(read2) 

    output:
    tuple val(meta), path("${meta}_sorted.bam"), emit: bam
    path "versions.yml"                        , emit: versions

    script:
    """
    bwa mem ${index} ${read1} ${read2} | samtools sort -o "${meta}_sorted.bam" -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	bwa: \$("bwa -version")
        samtools: \$(samtools version | head -n 1 | cut -d" " -f2)
    END_VERSIONS
    """
}

//below is alternate version of module for snps_phylo
process ALIGN_TO_REF_UBAM {

    label 'medium_task'
    label 'bwa'

    conda (params.enable_conda ? "bioconda::bwakit=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1':
        null }"

    input:
    path fastq
    tuple val(sampleId), file(reads)
    path ubam
    tuple val(meta), path(reference)
    tuple val(meta), path(reference_index)

    output:
    path "${fastq.simpleName}.bam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    bwa mem \
        -M \
        -t ${task.cpus} \
        -p ${reference} \
        ${fastq} > ${fastq.simpleName}.bam
    """
}


