process NANO_COMPARE_FASTQ {
    
    label 'medium_task'
    label 'nanocomp'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.12.0--py_0':
        null }"
    
    input:
    path reads
    
    output:
    path "compare_fastq", emit: fastq_graphs
    path("versions.yml"), emit: versions
    
    script:
    """
    NanoComp --fastq ${reads} \
        --outdir compare_fastq \
        -t ${task.cpus}

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(NanoComp -v | cut -d" " -f2)
    VERSIONS
    """
}