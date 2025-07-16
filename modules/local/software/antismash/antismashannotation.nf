process ANTISMASH_ANNOTATION {
    
    label 'antismash'
    label 'very_large_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::antismash=7.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'antismash/standalone:7.1.0': null }"
	
    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${meta}/*.gbk")    , emit: antismash_results

    script:
    """
    mkdir mpl
    export MPLCONFIGDIR=mpl
    antismash \\
        -c ${task.cpus} \\
        --output-dir ${meta} \\
        --genefinding-tool none \\
        --taxon fungi \\
        --tigrfam \\
        --cb-general \\
        --cb-knownclusters \\
        --pfam2go \\
        --rre \\
        --cc-mibig \\
        --clusterhmmer \\
        ${gbk}
    
    rm ${meta}/*.region*.gbk
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
