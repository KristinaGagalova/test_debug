process EXTRACT_KRAKEN_READS {

    label 'krakentools'
    label 'small_task'

    conda (params.enable_conda ? "bioconda::krakentools=1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0':
        null }"

    input:
    tuple val(meta), path(reads1), path(reads2), val(contam), path(kraken_report), path(kraken_out)

    output:
    tuple val(meta), path("${meta}_1_kraken2clean.fasta"), path("${meta}_2_kraken2clean.fasta"), emit:cleaned_reads 

    script:
    """
    python3 $projectDir/bin/extract_kraken_reads.py \\
    --exclude \\
    -s1 ${reads1} \\
    -s2 ${reads2} \\
    -r ${kraken_report} \\
    -k ${kraken_out} \\
    -o ${meta}_1_kraken2clean.fasta \
    -o2 ${meta}_2_kraken2clean.fasta \
    -t ${task.cpus}
    """
}