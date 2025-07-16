process FORMAT_READ_HEADERS {

    label 'small_task'

    // Input values
    input:
    tuple val(meta), path(reads) // list of fasta files [fasta1, fasta2] from channel collect()

    // Expected output
    output:
    tuple val(meta), path('formatted_r*.fastq.gz'), emit: formatted

    script:
    """
    if [[ ${reads[1]} == *.gz ]]; then
        gunzip -c ${reads[0]} > formatted_r1.fastq
        gunzip -c ${reads[1]} > formatted_r2.fastq
    else
        cp ${reads[0]} formatted_r1.fastq
        cp ${reads[1]} formatted_r2.fastq
    fi

    #sed -i 's:\\(.*\\).* \\(.*\\)\\( length=.*\\):\\1-\\2/1:g' formatted_r1.fastq
    #sed -i 's:\\(.*\\).* \\(.*\\)\\( length=.*\\):\\1-\\2/2:g' formatted_r2.fastq

    gzip formatted_r1.fastq
    gzip formatted_r2.fastq
    """
}
