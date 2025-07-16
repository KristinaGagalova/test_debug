process CONCAT_FASTA {

    tag "merging"
    label 'small_task'

    input:
    val(transcripts) // list of fasta files [fasta1, fasta2] from channel collect()

    output:
    tuple val("merged"), path('transcripts.merged.fasta'), emit: merged
    script:
    """
    cat ${transcripts.join(' ')} > transcripts.merged.fasta
    """
}

process CONCAT_ANNOTATIONS {

    tag 'merging'
    label 'small_task'

    input:
    tuple val(meta), path(funannotate_fasta), path(cqpm_fasta)

    output:
    tuple val(meta), path(funannotate_fasta)

    """
    cat ${funannotate_fasta} ${cqpm_fasta} > merged.fasta
    mv merged.fasta ${funannotate_fasta}
    """
}