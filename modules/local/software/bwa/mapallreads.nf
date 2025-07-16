process BWAMEM_MAPALL {

    label 'medium_task'
    label 'bwa'
    tag "Map reads ${sample_id} unpaired."

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    	'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_0'  :
        'quay.io/biocontainers/bwa:0.7.18--he4a0461_0'}"

    input:
    tuple val(sample_id), path(reads), path(index)     // reads from pair with sample id
    val(read_pair)                                     // either "R1" or "R2" from read pairs

    output:
    tuple val("${sample_id}"), path("${sample_id}_${read_pair}.sam"), emit: reads_mapped
    path("versions.yml")                                            , emit: versions

	
    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    bwa mem \\
	-t ${task.cpus} \\
	-a \$INDEX \\
	${reads} > ${sample_id}_${read_pair}.sam

    cat <<-VERSIONS > versions.yml
    "${task.process}":
         bwa: \$(bwa 2>&1 | grep \"Version\" | awk '{print \$2}')
    VERSIONS
    """
}