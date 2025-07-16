process EGGNOG_SETUP {

    label 'eggnog'
    label 'small_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0':
        null }"

    input:
    val eggnog_db

    output:
    path "versions.yml"       

    script:
    """
    if [ ! -d "${eggnog_db}" ]; then
        mkdir -p "${eggnog_db}"
    fi

    download_eggnog_data.py -y --data_dir ${eggnog_db}

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(emapper.py --version | grep "emapper" | sed -e "s/.*emapper-//" -e "s/ .*//")
    VERSIONS
    """
}
