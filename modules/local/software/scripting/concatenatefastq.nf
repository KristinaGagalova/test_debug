process CONCAT_FASTQ {

    tag "merging"
    label 'small_task'

    // Input values
    input:
    tuple val(meta), path(reads1), path(reads2) // list of fasta files [fasta1, fasta2] from channel collect()

    // Expected output
    output:
    tuple val(meta), path("${meta}_R1_combined.fastq.gz"), path("${meta}_R2_combined.fastq.gz"), emit: merged

    script:
    def r1_cmds = reads1.collect { it.toString().endsWith('.gz') ? "zcat $it" : "cat $it" }.join(' ')
    def r2_cmds = reads2.collect { it.toString().endsWith('.gz') ? "zcat $it" : "cat $it" }.join(' ')

    """
    ${r1_cmds} | gzip > ${meta}_R1_combined.fastq.gz
    ${r2_cmds} | gzip > ${meta}_R2_combined.fastq.gz
    """
}
