process GENERATE_TREE {

    label 'medium_task'
    label 'iqtree'

    conda (params.enable_conda ? "bioconda::iqtree=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.4--h21ec9f0_0':
        null }"

    input:
    path phylip

    output:
    path "*.treefile"

    """
    iqtree \
        -nt ${task.cpus}\
        -s ${phylip} \
        -st DNA \
        -m TEST \
        -bb 1000 \
        -alrt 1000
    """
}