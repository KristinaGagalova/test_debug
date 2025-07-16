# IMPORTANT: singularity cache dir needs to be specified in order to download image when running pante2 default
export NXF_SINGULARITY_CACHEDIR="./work"

nextflow run ./main.nf -profile pawsey_setonix,singularity \
		--tool pasa --genomes '../genomes/*fasta' \
		--transcripts '../transcripts/*fasta' \
		--outputdir testpasa -resume \
		--enable_conda false
