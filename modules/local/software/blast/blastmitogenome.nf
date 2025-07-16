process BLAST_MITOGENOME {

    label 'blast'
    label 'small_task'
    conda (params.enable_conda ? "bioconda::blast=2.16.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/blast:2.15.0--pl5321h6f7f691_1':
        null }"

    input:
    tuple val(meta), path(ncl_assembly), path(mito_assembly)

    output:
    tuple val(meta), path("${meta}_blast_hits_sorted.tsv"), emit:blast_table 

    script:
    """
    makeblastdb -in "${ncl_assembly}" -out "${meta}" -dbtype 'nucl' -hash_index

    blastn -query "${mito_assembly}" -task blastn -db "${meta}" -outfmt 6 -out "${meta}_blast_hits.tsv" -evalue 1e-10 -num_threads 16

    sort -n -r -k 12 "${meta}_blast_hits.tsv" > "${meta}_blast_hits_sorted.tsv"
    """
}