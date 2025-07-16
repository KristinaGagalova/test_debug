process STAR_MAPREADS {
    
    label 'star'
    label 'medium_task'
    tag "Star align reads for ${sample_id}"

    conda (params.enable_conda ? "bioconda::star=2.7.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0':
        'quay.io/biocontainers/star:2.7.9a--h9ee0642_0' }"
	
    output:
        tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"), emit: alignements
        tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Log.final.out")                , emit: reports
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_SJ.out.tab")                   , emit:splicejunctions

    input:
        tuple val(sample_id), path(reads)
		tuple val(sample_id), path(index_dir)
     
    script:
    """
    STAR --runThreadN ${task.cpus} \
		--runMode alignReads \
		--readFilesCommand zcat \
		--outSAMunmapped Within \
		--outFilterType BySJout \
		--outSAMattributes NH HI AS NM MD \
		--outSAMtype BAM SortedByCoordinate \
		--outTmpDir star_aligned/${sample_id}/_STARtmp \
		--outFileNamePrefix star_aligned/${sample_id}/${sample_id}_ \
		--genomeDir $index_dir \
		--readFilesIn ${reads[0].join(",")} ${reads[0].join(",")}
    """
}
