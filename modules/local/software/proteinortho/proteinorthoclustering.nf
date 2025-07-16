process PROTEIN_ORTHO_CLUSTER {

    label 'protein_ortho'
    label 'very_large_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::proteinortho:6.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
	"https://depot.galaxyproject.org/singularity/proteinortho:6.3.1--h70414c8_0":
        'quay.io/biocontainers/proteinortho:6.3.1--h70414c8_0'}"
		
    input:
    tuple path(prots), path(cds)

    output:
    path("myproject.proteinortho.tsv"), emit: protein_ortho

    script:
    """
    proteinortho -singles -selfblast -cpus=${task.cpus} ${prots}/*

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        proteinortho \$(proteinortho -version)
    VERSIONS
    """
}
