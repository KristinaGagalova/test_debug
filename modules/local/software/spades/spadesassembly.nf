process SPADES_ASSEMBLY {

    label 'spades'
    label 'very_large_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::spades=4.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:4.0.0--h5fb382e_1':
        null }"

    input: 
    tuple val(meta), path(unmerged_reads1), path(unmerged_reads2), path(merged_reads)

    output:
    tuple val(meta), path("${meta}_ncl_assembly.fasta"), emit: assembly
    tuple val(meta), path("${meta}"), emit: outdir

    """
    spades.py \\
        -1 ${unmerged_reads1} \\
        -2 ${unmerged_reads2} \\
        --merged ${merged_reads} \\
        --only-assembler \\
        --cov-cutoff auto \\
        -t ${task.cpus} \\
        -o ${meta}

    if [ -e "${meta}/scaffolds.fasta" ]; then
        cp ${meta}/scaffolds.fasta ${meta}_ncl_assembly.fasta
    else
        cp ${meta}/contigs.fasta ${meta}_ncl_assembly.fasta
    # Add your commands for when the file does not exist here
    fi
    """

}