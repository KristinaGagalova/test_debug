process FUNANNOTATE_SORT {
    label 'funannotate'
    label 'medium_task'
    tag "${meta}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
		
    input:
    tuple val(meta), path(fasta)
	
    output:
    tuple val(meta), path(fasta_clean), emit : fasta
    path("versions.yml")              , emit : versions

    script:
    fasta_clean = "${meta}.clean.sort.fa"
    """

    # Sort contigs using awk
    sortFasta.sh ${fasta} $fasta_clean

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS	
    """	
}
