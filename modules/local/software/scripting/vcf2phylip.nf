process VCF_TO_PHYLIP {

    label 'small_task'
    label 'python'

    conda (params.enable_conda ? "conda-forge::python=3.12.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/python:3.12':
        null }"

    input:
    path vcf

    output:
    path "${params.refname}.*.fasta"

    """
    python3 $projectDir/bin/vcf2phylip.py \
	    -i ${vcf} \
	    --fasta \
	    --nexus
    """
}