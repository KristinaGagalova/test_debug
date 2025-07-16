process CQ_TRANSCRIPTOME_PREDICT {

    label 'funannotate'
    label 'biggish_task'
    tag "${meta}"
    errorStrategy 'ignore'

    conda (params.enable_conda ? "bioconda::funannotate=1.8.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nextgenusfs/funannotate:v1.8.15':
        'https://depot.galaxyproject.org/singularity/funannotate:1.8.15--pyhdfd78af_2' }"
			
    input:
    tuple val(meta), path(fasta_genome), path(gff)

    output:
    tuple val(meta), path(fasta_genome), path(gff), path("${meta}_PredictedPass.gff3")

    """
    CodingQuarry -f ${fasta_genome} -t ${gff} -p ${task.cpus}
    mv out/PredictedPass.gff3 ${meta}_PredictedPass.gff3

    rm -rf ParameterFiles
    """
}

//---------------------------------------------------
//-------- Processes for processing for Pathogen mode
//---------------------------------------------------

process CQ_EXTRACT_CDS {
    label 'coding_quarry'
    label 'small_task'
    tag "${meta}"

    container "docker://pavmis/cqpm:1.1"

    input:
    tuple val(meta), path(fasta_genome), path(gff)

    output:
    tuple val(meta), path("PM_${meta}_PredictedPass.fa"), emit: predicted_cds
    
    script:
    """	
    CodingQuarry -f ${fasta_genome} -a ${gff} -1 -g none
    mv out/Predicted_CDS.fa PM_${meta}_PredictedPass.fa
    """
}

process CQ_CDS_FILTER_PROCESS {
    label 'coding_quarry'
    label 'small_task'
    tag "${meta}"

    container "docker://pavmis/cqpm:1.1"
    
    input:    
    tuple val(meta), path(extracted_cds)

    output:
    tuple val(meta), path("CQ_${meta}_Proteins.fa-*"), emit: split_files

    script:
    """
    python \$QUARRY_PATH/scripts/fastaTranslate.py ${extracted_cds} | sed 's/*\$//g' > CQ_Proteins.fa
    python \$QUARRY_PATH/scripts/gene_errors_Xs.py CQ_Proteins.fa CQPMtmp.fa
    mv CQPMtmp.fa CQ_${meta}_Proteins.fa

    # split the protein file up into smaller chunks
    python \$QUARRY_PATH/scripts/split_fasta.py CQ_${meta}_Proteins.fa 200
    """
}

process CQ_PATHOGENMODE {
    errorStrategy 'ignore'
    tag "${meta}"
    label 'coding_quarry'
    label 'biggish_task'

    container "docker://pavmis/cqpm:1.1"

    input:
    tuple val(meta), path(gff), path(assembly), path(secretome)

    output:
    tuple val(meta), path("${meta}_CQPM.gff3"), emit: output
    tuple val(meta), path("${meta}_CQPM.cds.fa"),    emit: cds
    tuple val(meta), path("${meta}_CQPM.prot.fa"),   emit: protseq

    script:
    """
    CodingQuarry -a ${assembly} -f ${gff} -p ${task.cpus} -g ${secretome} -h

    mv out/PGN_predictedPass.gff3 ${meta}_CQPM.gff3
    mv out/PGN_predicted_CDS.fa ${meta}_CQPM.cds.fa

    python \$QUARRY_PATH/scripts/fastaTranslate.py ${meta}_CQPM.cds.fa | sed 's/*\$//g' > ${meta}_CQPM.prot.fa
    """
}
