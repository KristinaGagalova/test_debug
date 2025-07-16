process FUNANNOTATE_TRAIN {
    label 'funannotate'
    label 'big_task_highmem'
    scratch true

    tag "${meta}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
		
    input:
    tuple val(meta), path(fasta_genome), path(rnareads) // input genome [meta, fasta, rna-seq]
    val(species)


    output:
    tuple val(meta), path("${meta}_training"),           emit: training
    path("${meta}_transcripts.fasta"),          emit: transcripts
    tuple val(meta), path("${meta}_rna.bam"),   emit: bam

    script:
    """
    funannotate train \\
        -i      ${fasta_genome} \\
        -o      out \\
        -l      ${rnareads[0]} \\
        -r      ${rnareads[1]} \\
        --cpus  ${task.cpus} \\
        --species ${species}
    
    mv out/training/trinity.fasta ${meta}_transcripts.fasta
    mv out/training/hisat2.coordSorted.bam ${meta}_rna.bam
    mv out ${meta}_training


    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}
