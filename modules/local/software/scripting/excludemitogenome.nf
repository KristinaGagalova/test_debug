process EXCLUDE_MITOGENOME {

    label 'python'
    label 'small_task'

    conda (params.enable_conda ? "conda-forge::python=3.12.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/python:3.12':
        null }"

    input:
    tuple val(meta), path(blast_hits), path(assembly), path(index)

    output:
    tuple val(meta), path("${meta}_ncl_assembly_filtered.fasta"), emit:excluded_fasta

    script:
    """
    calculate_length.py ${blast_hits}

    exclude_overlapped.py ${meta}_lengths.tsv ${index}

    f="${assembly}"
    exclude="${meta}_exclude.txt"
    output="${meta}_ncl_assembly_filtered.fasta"
    awk '{ if ((NR>1)&&(\$0~/^>/)) { printf("\\n%s", \$0); } else if (NR==1) { printf("%s", \$0); } else { printf("\\t%s", \$0); } }' "\$f" | grep -vFf "\$exclude" - | tr "\\t" "\\n" > "\$output"
    """
}