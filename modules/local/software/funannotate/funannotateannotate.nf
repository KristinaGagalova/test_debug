process FUNANNOTATE_ANNOTATE {
    label 'funannotate'
    label 'big_task'
    tag "${meta}"

    conda (params.enable_conda ? "bioconda::funannotate=1.8.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
		
    input:
    tuple val(meta), path(assembly), path(gff), path(eggnog_results), path(antismash_results)
    val(species)

    output:
    path("${meta}")                                                             , emit: output
    tuple val(meta), path("${meta}/annotate_results/*.gff3")                    , emit: gff

    script:
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
    def funannotate_db = params.funannotate_db == null ? "${baseDir}/work/funannotate_db" : "${baseDir}/work/funannotate_db"
    dev args_header_len = params.header_length ? "--header_length '${params.header_length}'" : ""
    """
    export FUNANNOTATE_DB=$funannotate_db

    funannotate annotate \\
        --gff ${gff} \\
        --fasta ${assembly} \\
        --out functional_annotations \\
        --species ${species} \\
        --cpus ${task.cpus} \\
        --eggnog ${eggnog_results} \\
        --antismash ${antismash_results} \\
        ${args_buscodb} \\
	${ args_header_len}
    
    mv functional_annotations ${meta}

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}
