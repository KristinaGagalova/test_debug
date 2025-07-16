process ASSEMBLY_STATS {

    label 'small_task'

    // Input values
    input:
    tuple val(meta), path(ncl_assembly), path(mito_assembly)  

    output:
    tuple val(meta), path("${meta}_spades.txt"), emit: spades_stats
    tuple val(meta), path("${meta}_mitoz.txt"), emit: mitoz_stats

    script:
    """
    n50calc_singleline.pl ${ncl_assembly} 200 >> ${meta}_spades.txt
    n50calc_singleline.pl ${mito_assembly} 200 >> ${meta}_mitoz.txt
    """
}
