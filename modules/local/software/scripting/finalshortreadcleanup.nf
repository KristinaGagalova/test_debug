process FINAL_SHORT_READ_CLEANUP {
    
    label 'bioawk'
    label 'small_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::bioawk=1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6':
        null }"

    input:
    tuple val(meta), path(ncl_assembly), path(mito_assembly)

    output:
    tuple val(meta), path("${meta}_final_assembly.fasta"), emit:final_assembly

    script:
    """
    bioawk -c fastx -v prefix="${meta}" '{if(length(\$seq)>=200) print ">" prefix "_scaffold-" (++i)" "length(\$seq)"\\n"\$seq }' < "${ncl_assembly}" > "${meta}_ncl_assembly_filtered_sorted.fasta"

    bioawk -c fastx -v prefix="${meta}" '{ print ">" prefix "_mito-" (++i)" "length(\$seq)"\\n"\$seq }' < "${mito_assembly}" > "${meta}_mito_assembly_sorted.fasta"

    cat "${meta}_ncl_assembly_filtered_sorted.fasta" "${meta}_mito_assembly_sorted.fasta" > "${meta}_final_assembly.fasta"
    """
}

process FINAL_SHORT_READ_CLEANUP_NO_MITOZ {

    label 'bioawk'
    label 'small_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::bioawk=1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/bioawk:1.0--he4a0461_12':
        null }"

    input:
    tuple val(meta), path(ncl_assembly)

    output:
    tuple val(meta), path("${meta}_final_assembly.fasta"), emit:final_assembly

    script:
    """
    bioawk -c fastx -v prefix="${meta}" '{if(length(\$seq)>=200) print ">" prefix "_scaffold-" (++i)" "length(\$seq)"\\n"\$seq }' < "${ncl_assembly}" > "${meta}_ncl_assembly_filtered_sorted.fasta"

    cat "${meta}_ncl_assembly_filtered_sorted.fasta" > "${meta}_final_assembly.fasta"
    """
}