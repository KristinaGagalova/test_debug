process HELPER_PASA2TRAINING {

    tag "${meta}"
    label 'pasa'
    label 'process_low'
    
    container 'pasapipeline/pasapipeline:2.5.2'
    
    input:
    tuple val(meta), path(gff)
    val(nmodels)

    output:
    tuple val(meta), path(training_gff), emit: gff
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    training_gff = prefix + ".pasa_training.gff3"
    """
    pasa_select_training_models.pl --nmodels $nmodels --infile $gff >> $training_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        helper: 1.0
    END_VERSIONS
    """
}
