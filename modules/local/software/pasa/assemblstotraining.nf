process PASA_ASMBLSTOTRAINING {

    tag "${meta}"
    label 'pasa'
    label 'process_medium'
    
    container 'pasapipeline/pasapipeline:2.5.2'
    //container 'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2'

    input:
    tuple val(meta), path(fasta), path(gff)

    output:
    tuple val(meta), path('*.genome.gff3')     , emit: gff
    tuple val(meta), path('*.transdecoder.pep'), emit: fasta
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    \$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi \
       --pasa_transcripts_fasta $fasta \
       --pasa_transcripts_gff3 $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
		pasa: 2.5.2
    END_VERSIONS
    """
}
