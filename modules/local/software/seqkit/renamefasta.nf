process SEQKIT_LISTRENAME {

    label 'small_task'
    label 'seqkit'
    tag "Rename final assembly ${sample_id}."

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0'  :
        'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'}"

    input:
    tuple val(sample_id), path(fasta_assembly)

    output:
    tuple val("${sample_id}"), path("${sample_id}.tab"), emit: list_scafs
	path("versions.yml")                               , emit: versions

    script:
    """
    seqkit sort -lr ${fasta_assembly} |\
	grep ">" | sed 's/>//' |\
	awk 'BEGIN {FS=OFS="\t"} {print \$1,NR}' |\
	awk 'BEGIN {FS=OFS="\t"} {\$2 = sprintf("%02d", \$2)}1' |\
	awk -v s=${sample_id} 'BEGIN {FS=OFS="\t"} {print \$1,s"_scaf"\$2}' > ${sample_id}.tab

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    VERSIONS
    """

}


process SEQKIT_RENAMEGEN {
 
    label 'small_task'
    label 'seqkit'
    tag "Rename final assembly ${sample_id}."
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0'  :
        'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'}"

    input:
	tuple val(sample_id), path(list_nams), path(fasta_assembly)

    output:
	tuple val(sample_id), path("${sample_id}_renamed.fasta"), emit: renamed_fasta
	path("versions.yml")                                    , emit: versions
    
    script:
    """
    seqkit replace \
	-p '(.+)\$' \
	-r '{kv}'\
	-k ${list_nams} \
	${fasta_assembly} 1> ${sample_id}_renamed.fasta
    
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    VERSIONS
    """
}


process SEQKIT_MITOTAG {

    label 'small_task'
    label 'seqkit'
    tag "Tag mito genome in assembly ${sample_id}."

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0'  :
        'quay.io/biocontainers/seqkit:2.8.1--h9ee0642_0'}"

    input:
    tuple val(sample_id), path(list_mito), path(fasta_assembly)

    output:
    tuple val(sample_id), path("${sample_id}_final.fasta"), emit: mito_fasta 
	path("versions.yml")                                  , emit: versions

    script:
    """
    seqkit replace \
	--keep-key \
	-p '(.+)\$' \
	-r '{kv}' \
	-k ${list_mito} \
	${fasta_assembly} 1> ${sample_id}_final.fasta
    
    cat <<-VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    VERSIONS
    """
}