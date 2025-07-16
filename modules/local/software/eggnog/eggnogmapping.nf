process EGGNOG_MAPPING {

    label 'eggnog'
    label 'medium_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    	'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_2':
        'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0'}"
	
    input:
    tuple val(meta), path(protseq)
    val(eggnog_db)

    output:
    tuple val(meta), path("*.emapper.annotations")  , emit: eggnog_results

    script:
    """
    emapper.py  -i ${protseq} \
                -o ${meta} \
                --cpu ${task.cpus} \
                --data_dir ${eggnog_db}
    
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(emapper.py --version | grep "emapper" | sed -e "s/.*emapper-//" -e "s/ .*//")
    VERSIONS
    """
}
