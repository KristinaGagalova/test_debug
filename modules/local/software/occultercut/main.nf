/*
 * OcculterCut
 * doi: 10.1093/gbe/evw121
 */
process OCCULTERCUT {
    label "small_task"
    time "4h"
    tag "${meta}"

    conda "bioconda::occultercut=1.1"

	//not working - TODO fix
    //conda "${software.occultercut.conda_channel}::occultercut=${software.occultercut.version}"
   
    publishDir "${params.outdir}/${meta}/noncoding", mode: 'copy'

    input:
    tuple val(meta), file(fasta)

    output:
    tuple val(meta), file("${meta}_occultercut_regions.gff3"), emit: occulterCutRegions
    tuple val(meta), file("${meta}_occultercut_grouped_regions.gff3"), optional: true, emit: occulterCutGroupedRegions
    path "${meta}_occultercut.svg", optional: true
    path "${meta}_occultercut_composition_gc.txt"
    path "${meta}_occultercut_my_genome.txt"
    path "${meta}_occultercut_nuc_frequencies.R*", optional: true
   
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    OcculterCut -f "${fasta}"

    if [ -e plot.plt ]
    then
      sed -i '1i set terminal svg size 900,600 enhanced font "Helvetica,20"' plot.plt
      sed -i '1a set output "plot.svg"' plot.plt
      gnuplot plot.plt
      mv plot.svg "${meta}_occultercut.svg"
    fi

    mv compositionGC.txt "${meta}_occultercut_composition_gc.txt"
    mv regions.gff3 "${meta}_occultercut_regions.gff3"
    mv myGenome.txt "${meta}_occultercut_my_genome.txt"

    if [ -e groupedRegions.gff3 ]
    then
      mv groupedRegions.gff3 "${meta}_occultercut_grouped_regions.gff3"
    fi

    find . -name "nuc_frequences.R*" -printf '%f\\0' \
    | xargs -I {} -0 -- mv '{}' "${meta}_occultercut_{}"
    """
}