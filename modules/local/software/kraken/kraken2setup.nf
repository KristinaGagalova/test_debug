process KRAKEN2_SETUP {

    label 'kraken'
    label 'very_large_task'

    conda (params.enable_conda ? "bioconda::kraken2=2.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_1':
        null }"

    input:
    val(db)

    output:
    val(db_complete)

    script:
    """
    kraken2-build --standard --db ${db} --threads ${task.cpus}
    """
}