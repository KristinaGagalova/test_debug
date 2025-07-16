process FUNANNOTATE_SETUP {
    label 'funannotate'
    label 'small_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::funannotate=1.8.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
		
    input:
    val buscodb
    val funannotate_db

    output:
    path "versions.yml"         , emit: versions

    script:
    """
    funannotate setup --wget -b $buscodb -d $funannotate_db
	
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}
