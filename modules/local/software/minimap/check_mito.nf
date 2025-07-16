process MITO_ALIGN {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                    'https://depot.galaxyproject.org/singularity/pomoxis%3A0.3.15--pyhdfd78af_0':
                    'quay.io/biocontainers/pomoxis:0.2.2--py_0' }"

    label 'small_task'
    label 'pomoxis'
    tag "Get mito genome ${sample_id}"

    input:
	tuple val(sample_id), path(fasta_assembly)
	file(mitodb)

    output:
    tuple val(sample_id), path("${sample_id}_MT.tab"), emit: list_mito
	path("versions.yml")                             , emit: versions

    script:
    """
    minimap2 -ax asm5 ${mitodb} ${fasta_assembly} |\
	samtools view -F 4 | awk '{print \$1,\$1,"MT"}' | sed 's/ /\t/' | sort | uniq 1> ${sample_id}_MT.tab

    cat <<-VERSIONS > versions.yml
    "${task.process}":
	 minimap2: \$(minimap2 --version 2>&1)
    VERSIONS    
    """
}