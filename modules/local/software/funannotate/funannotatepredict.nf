process FUNANNOTATE_GENE_PREDICTION {
    label 'funannotate'
    maxForks 20
    label 'big_task'
    tag "${meta}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
		
    input:
    val(versions)
    tuple val(meta), path(fasta_genome) // input genome [(meta), (fasta)]
    val(protein)
    tuple val(meta_trans), path(fasta_trans) // input transcriptome [(meta), (fasta)]
    val(species)
    val(rnabam)
    val(stringtie)

    output:
    path("${meta}")                                                 , emit: output
    tuple val(meta), path("${meta}.proteins.fa")                    , emit: protseq
    tuple val(meta), path("${meta}.cds.fa")                         , emit: cds
    tuple val(meta), path("${meta}/predict_results/*.gff3")         , emit: gff
    tuple val(meta), path("${meta}/predict_results/*.gbk")          , emit: gbk
    path "versions.yml"                                             , emit: versions

    script:
    def funannotate_db = params.funannotate_db == null ? "${baseDir}/work/funannotate_db" : "${baseDir}/work/funannotate_db"
    def genmark_path = params.gmark_db
    // def augustus_config = params.augustus_config
    // def evm_home = params.evm_home
    def args_prot = protein != "EMPTY" ? "--protein_evidence '${protein}'" : ""
    def args_trans = fasta_trans != "EMPTY" ? "--transcript_evidence '${fasta_trans}'" : ""
    def args_rna = rnabam != "EMPTY" ? "--rna_bam '${rnabam}'" : ""
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
    def args_buscoseed = params.buscoseed == null ? "" : "--busco_seed_species '${params.buscoseed}'"
    def args_ploidy = params.ploidy ? "--ploidy '${params.ploidy}'" : ""
    """
    export FUNANNOTATE_DB=$funannotate_db
    export GENEMARK_PATH=$genmark_path
    augustus_species=\$(echo ${species} | sed 's/ /_/g')

    funannotate predict -i ${fasta_genome} \\
    	-o ${meta} \\
    	-s "${species}" \\
    	--augustus_species \$augustus_species \\
    	--cpus ${task.cpus} \\
    	--organism "$params.organism" \\
    	${args_buscodb} \\
     	${args_ploidy} \\
        --min_training_models 100\\
     	${args_buscoseed} ${args_prot} ${args_trans} ${args_rna}

    cp ${meta}/predict_results/*.proteins.fa ${meta}.proteins.fa
    cp ${meta}/predict_results/*.cds-transcripts.fa ${meta}.cds.fa

    rm -rf ${meta}/predict_misc

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}


process FUNANNOTATEPASA_GENE_PREDICTION_CQ {
    label 'funannotate'
    maxForks 20
    label 'big_task'
    tag "${meta}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"

    input:
    val(versions)
    tuple val(meta), path(fasta_genome), path(pasa_gff), path(cq_gff) // input genome [(meta, fasta, pasa_gff),...]
    val(protein)
    tuple val(meta_trans), path(fasta_trans) // input transcriptome [(meta), (fasta)]
    val(species)
    //val(stringtie)

    output:
    path("${meta}")                                                 , emit: output
    tuple val(meta), path("${meta}.proteins.fa")                    , emit: protseq
    tuple val(meta), path("${meta}.cds.fa")                         , emit: cds
    tuple val(meta), path("${meta}/predict_results/*.gff3")         , emit: gff
    tuple val(meta), path("${meta}/predict_results/*.gbk")          , emit: gbk
    path "versions.yml"                                             , emit: versions

    script:
    def funannotate_db = params.funannotate_db == null ? "${baseDir}/work/funannotate_db" : "${baseDir}/work/funannotate_db"
    def genmark_path = params.gmark_db
    // def augustus_config = params.augustus_config
    // def evm_home = params.evm_home 
    def args_prot = protein != "EMPTY" ? "--protein_evidence '${protein}'" : ""
    def args_trans = fasta_trans != "EMPTY" ? "--transcript_evidence '${fasta_trans}'" : ""
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
    def args_buscoseed = params.buscoseed == null ? "" : "--busco_seed_species '${params.buscoseed}'"
    def args_ploidy = params.ploidy ? "--ploidy '${params.ploidy}'" : ""
    def args_cq = cq_gff != "EMPTY" ? "--other_gff '${cq_gff}:2'" : ""
    
    """
    export FUNANNOTATE_DB=$funannotate_db
    export GENEMARK_PATH=$genmark_path
    augustus_species=\$(echo ${species} | sed 's/ /_/g')

    funannotate predict -i ${fasta_genome} \\
        -o ${meta} \\
        -s "${species}" \\
        --augustus_species \$augustus_species \\
        --cpus ${task.cpus} \\
        --organism "$params.organism" \\
        ${args_buscodb} \\
        ${args_ploidy} \\
	    --pasa_gff ${pasa_gff} \\
        --min_training_models 100 \\
        ${args_buscoseed} ${args_prot} ${args_trans} \\
        ${args_cq}\\
        --weights codingquarry:0
    
    cp ${meta}/predict_results/*.proteins.fa ${meta}.proteins.fa
    cp ${meta}/predict_results/*.cds-transcripts.fa ${meta}.cds.fa

    rm -rf ${meta}/predict_misc

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}

process FUNANNOTATEPASA_GENE_PREDICTION_BAM {
    label 'funannotate'
    maxForks 20
    label 'big_task'
    tag "${meta}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"

    input:
    val(versions)
    tuple val(meta), path(fasta_genome), path(rnabam), path(pasa_gff) // input genome [meta, fasta, pasa_gff,...]
    val(protein)
    tuple val(meta_trans), path(fasta_trans) // input transcriptome [(meta), (fasta)]
    val(species)
    //val(stringtie)

    output:
    path("${meta}")                                                 , emit: output
    tuple val(meta), path("${meta}.proteins.fa")                    , emit: protseq
    tuple val(meta), path("${meta}.cds.fa")                         , emit: cds
    tuple val(meta), path("${meta}/predict_results/*.gff3")         , emit: gff
    tuple val(meta), path("${meta}/predict_results/*.gbk")          , emit: gbk
    path "versions.yml"                                             , emit: versions

    script:
    def funannotate_db = params.funannotate_db == null ? "${baseDir}/work/funannotate_db" : "${baseDir}/work/funannotate_db"
    def genmark_path = params.gmark_db
    // def augustus_config = params.augustus_config
    // def evm_home = params.evm_home
    def args_prot = protein != "EMPTY" ? "--protein_evidence '${protein}'" : ""
    def args_trans = fasta_trans != "EMPTY" ? "--transcript_evidence '${fasta_trans}'" : ""
    def args_rna = rnabam != "" && rnabam != "EMPTY" ? "--rna_bam '${rnabam}'" : ""
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
    def args_buscoseed = params.buscoseed == null ? "" : "--busco_seed_species '${params.buscoseed}'"
    def args_ploidy = params.ploidy ? "--ploidy '${params.ploidy}'" : ""

    """
    export FUNANNOTATE_DB=$funannotate_db
    export GENEMARK_PATH=$genmark_path
    augustus_species=\$(echo ${species} | sed 's/ /_/g')

    funannotate predict -i ${fasta_genome} \\
        -o ${meta} \\
        -s "${species}" \\
        --augustus_species \$augustus_species \\
        --cpus ${task.cpus} \\
        --organism "$params.organism" \\
        ${args_buscodb} \\
        ${args_ploidy} \\
        --pasa_gff ${pasa_gff} \\
        --min_training_models 100 \\
        ${args_buscoseed} ${args_prot} ${args_trans} ${args_rna}\\
        --weights codingquarry:0
    
    cp ${meta}/predict_results/*.proteins.fa ${meta}.proteins.fa
    cp ${meta}/predict_results/*.cds-transcripts.fa ${meta}.cds.fa

    rm -rf ${meta}/predict_misc

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}


process FUNANNOTATE_GENE_PREDICTION_TRAINED {
    label 'funannotate'
    label 'big_task'
    tag "${meta}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"

    input:
    val(versions)
    tuple val(meta), path(fasta_genome), path(training_data), path(rnabam)
    val(protein)
    val(species)
 

    output:
    path("${meta}")                                                 , emit: output
    tuple val(meta), path("${meta}.proteins.fa")                    , emit: protseq
    tuple val(meta), path("${meta}.cds.fa")                         , emit: cds
    tuple val(meta), path("${meta}/predict_results/*.gff3")         , emit: gff
    tuple val(meta), path("${meta}/predict_results/*.gbk")          , emit: gbk
    path "versions.yml"                                             , emit: versions

    script:
    def funannotate_db = params.funannotate_db == null ? "${baseDir}/work/funannotate_db" : "${baseDir}/work/funannotate_db"
    def genmark_path = params.gmark_db
    // def augustus_config = params.augustus_config
    // def evm_home = params.evm_home
    def args_prot = protein != "EMPTY" ? "--protein_evidence '${protein}'" : ""
    def args_rna = rnabam != "" && rnabam != "EMPTY" ? "--rna_bam '${rnabam}'" : ""
    def args_buscodb = params.buscodb ? "--busco_db '${params.buscodb}'" : ""
    def args_buscoseed = params.buscoseed == null ? "" : "--busco_seed_species '${params.buscoseed}'"
    def args_ploidy = params.ploidy ? "--ploidy '${params.ploidy}'" : ""
    """
    export FUNANNOTATE_DB=$funannotate_db
    export GENEMARK_PATH=$genmark_path
    augustus_species=\$(echo ${species} | sed 's/ /_/g')

    funannotate predict -i ${fasta_genome} \\
        -o ${meta} \\
        -s "${species}" \\
        --cpus ${task.cpus} \\
        ${args_buscodb} \\
        ${args_ploidy} \\
        ${args_buscoseed} \\
        ${args_prot} \\
        ${args_rna} \\
        --augustus_species \$augustus_species

    cp ${meta}/predict_results/*.proteins.fa ${meta}.proteins.fa
    cp ${meta}/predict_results/*.cds-transcripts.fa ${meta}.cds.fa

    rm -rf ${meta}/predict_misc

    cat <<-VERSIONS > versions.yml
    "${task.process}":
        funannotate: \$(funannotate -version | sed -e "s/funannotate v//g")
    VERSIONS
    """
}
