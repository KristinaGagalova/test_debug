process getOcculterCutRegionFrequencies {
    label "small_task"
    time "2h"

    container 'git+https://github.com//darcyabjones/gffpal.git@454424857b7cc914d116c708ab73da095fa07fb2'	
	//does this work? need to confirm/fix

    tag "${meta}"

    input:
    tuple val(meta), file("in.gff"), file("in.fasta")

    output:
    tuple val(meta),
        val("occultercut_frequencies"),
        file("out.gff"), emit: occulterCutRegionFrequencies

    script:
    """
    tidy_occultercut_regions.py \
      --source "OcculterCut" \
      --type "region" \
      -o out.gff \
      in.gff \
      in.fasta
    """
}

process getOcculterCutGroupedRegionFrequencies {
    label "small_task"
    time "2h"

    tag "${meta}"

    input:
    tuple val(meta), file("in.gff"), file("in.fasta")

    output:
    tuple val(meta),
	val("occultercut_grouped_frequencies"),	
	file("out.gff"), emit: getOcculterCutGroupedRegionFrequencies

    script:
    """
    awk '
      BEGIN {OFS="\\t"; i=1}
      {print \$1, \$2, "region", \$4, \$5, \$6, \$7, \$8, "Name="\$3}
      ' in.gff \
    | tidy_occultercut_regions.py \
        --source "OcculterCut" \
        --type "region" \
        -o out.gff \
        - \
        in.fasta
    """
}

process tidyOcculterCutGFFs {
    label "small_task"
    time "1h"
    
    conda "bioconda::genometools-genometools=1.6.5"    
    
	//not working - TODO fix
    //conda "${software.genometools.conda_channel}::genometools-genometools=${software.genometools.version}"

    tag "${meta}"

    publishDir "${params.outdir}/${meta}/noncoding", mode: 'copy'

    input:
    tuple val(meta),
        val(suffix),
        file("in.gff")

    output:
    path "${meta}_${suffix}.gff3", emit: occulterCutGff

    script:
    """
    gt gff3 \
      -tidy \
      -sort \
      in.gff \
    > "${meta}_${suffix}.gff3"
    """
}