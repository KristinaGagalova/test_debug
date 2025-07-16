process FUNANNOTATE_CLEAN {
    label 'seqkit'
    label 'medium_task'
    tag "${meta}"
	
    conda (params.enable_conda ? "bioconda::seqkit=2.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0' }"
		
    input:
    tuple val(meta), path(fasta)
	
    output:
    tuple val(meta), path(fasta_clean), emit: fasta
    path "versions.yml"                   , emit: versions                     

    script:
    fasta_clean = fasta.getBaseName() + "_clean.fasta"
    min_contig = params.min_contig_size == null ? 500 : "${params.min_contig_size}"
    //removed command
    //funannotate clean -i $fasta -o $fasta_clean -m $min_contig
    //removed version:
    //funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    """
    seqkit seq -m $min_contig $fasta > $fasta_clean 

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version)
    VERSIONS
    """
}
