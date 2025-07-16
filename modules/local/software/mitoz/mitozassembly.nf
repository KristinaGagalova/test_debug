process MITOZ_ASSEMBLY {

    label 'mitoz'
    label 'large_task'
    tag "${meta}"
    errorStrategy 'ignore'

    conda (params.enable_conda ? "bioconda::mitoz=3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/mitoz:3.6--pyhdfd78af_1':
        null }"

    input: 
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("megahit/mitochondrion_final.contigs.fa.reformatted.fa")

    """
    extension="${reads1.Extension}"
    fastq_length=150
    if [ \$extension == "gz" ]; then
        fastq_length="\$(zcat ${reads1} | head -n 42 | sed -n '42p'  | wc -c)"
    else
        fastq_length="\$(sed -n '42p' ${reads1} | wc -c)"
    fi

    mitoz assemble \\
        --genetic_code 4 \\
        --clade Chordata \\
        --insert_size 250 \\
        --thread_number ${task.cpus} \\
        --fq1 ${reads1} \\
        --fq2 ${reads2} \\
        --fastq_read_length "\$fastq_length" \\
        --outprefix ${meta} \\
        --requiring_taxa Chordata \\
        --filter_by_taxa \\
        1>${meta}.log 2>${meta}.err
    """

}

process MITOZ_ASSEMBLY_SINGLE_END {

    label 'mitoz'
    label 'large_task'
    tag "${meta}"
    errorStrategy 'ignore'

    conda (params.enable_conda ? "bioconda::mitoz=3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/mitoz:3.6--pyhdfd78af_1':
        null }"

    input: 
    tuple val(meta), path(reads1)

    output:
    tuple val(meta), path("megahit/mitochondrion_final.contigs.fa.reformatted.fa")

    """
    extension="${reads1.Extension}"
    fastq_length=150
    if [ \$extension == "gz" ]; then
        fastq_length="\$(zcat ${reads1} | head -n 42 | sed -n '42p'  | wc -c)"
    else
        fastq_length="\$(sed -n '42p' ${reads1} | wc -c)"
    fi

    mitoz assemble \\
        --genetic_code 4 \\
        --clade Chordata \\
        --insert_size 250 \\
        --thread_number ${task.cpus} \\
        --fq1 ${reads1} \\
        --fastq_read_length "\$fastq_length" \\
        --outprefix ${meta} \\
        --requiring_taxa Chordata \\
        --filter_by_taxa \\
        1>${meta}.log 2>${meta}.err
    """

}