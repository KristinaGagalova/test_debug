process PREDICT_REPRESENTATIVES {
    label "metaeuk"
    label "medium_task"

    conda (params.enable_conda ? "bioconda::metaeuk:7.bba0d80--pl5321hd6d6fdc_2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    	'https://depot.galaxyproject.org/singularity/metaeuk:7.bba0d80--pl5321hd6d6fdc_2': 
    	'quay.io/biocontainers/metaeuk:7.bba0d80--pl5321hd6d6fdc_2'}"

    input:
    tuple val(meta), path(assemblies)
    path(protein_reps)

    output:
    path("${meta}_metaeuk.*"), emit: metaeuk_output_all
    path("${meta}_metaeuk.codon.fas"), emit: metaeuk_output_codon_fas
    path("${meta}_metaeuk.fas"), emit: metaeuk_output_fas
    path("${meta}_metaeuk.headersMap.tsv"), emit: metaeuk_output_tsv
    path("${meta}_metaeuk.gff"), emit: metaeuk_output_gff
		 
    script:
    """
    metaeuk easy-predict ${assemblies} ${protein_reps} ${meta}_metaeuk /tmp/${meta} --protein 1 --max-intron 500 --min-exon-aa 5 --max-overlap 5 --write-frag-coords 1
    """
}
