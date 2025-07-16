process KRAKEN2 {

    label 'kraken'
    label 'largeHighMem_task'

    conda (params.enable_conda ? "bioconda::kraken2=2.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1':
        null }"

    input:
    tuple val(meta), path(reads1), path(reads2), val(contam)
    val(kraken_db)

    output:
    tuple val(meta), path("${meta}_kraken2.report"), path("${meta}_kraken2.out"), emit:kraken_output

    script:
    """
    kraken2 \\
    --db ${kraken_db} \\
    --paired ${reads1} ${reads2} \\
    --threads ${task.cpus} \\
    --report ${meta}_kraken2.report \\
    --output ${meta}_kraken2.out
    """
}