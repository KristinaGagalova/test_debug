process STAR_INDEX {

    label 'star'
    label 'small_task'
    tag { "star: index" }

    conda (params.enable_conda ? "bioconda::star=2.7.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0':
        'quay.io/biocontainers/star:2.7.9a--h9ee0642_0' }"
		
    input:
    tuple val(meta), path(genome)
    
    output:
    tuple val(meta), path("${meta}_starindex"), emit: star_index
    
    script:
    """
    mkdir ${meta}_starIndex

    STAR --runMode genomeGenerate \
		--genomeSAindexNbases 11 \
		--genomeDir "${meta}_starindex" \
		--genomeFastaFiles ${genome} \
		--runThreadN ${task.cpus}
    """
}
