#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64GB #requirements may be high in case user downloads large images
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/23.10.0
module load singularity/4.1.0

# IMPORTANT: singularity cache dir needs to be specified in order to download image when running pante2 default
export NXF_SINGULARITY_CACHEDIR="./work"
mkdir -p singularity_tmpdir
export SINGULARITY_TMPDIR=$PWD/singularity_tmpdir

nextflow run ./main.nf -profile pawsey_setonix,singularity \
		-resume \
		-with-dag flowchart.png \
		--tool annotate-pangenome \
		--genomes '/scratch/y95/kgagalova/mp-cleanup/genomes_soft/*.fa' \
		--transcripts '/scratch/y95/kgagalova/mp-cleanup/transcripts/*fasta' \
		--rnareads '/scratch/y95/kgagalova/mp-cleanup/rnaseq/reads-run/*.{1,2}.fastq' \
		--outputdir results \
		--enable_conda false \
		--gmark_db --gmark_db "/scratch/y95/kgagalova/src/genemark/gmes_linux_64_4/" \
		--signalp "/scratch/y95/kgagalova/mp-cleanup/test/signalp-4.1g-work/signalp-4.1g-install/"
		# --augustus_config "/scratch/y95/kgagalova/software/Augustus/config" \
		# --evm_home "/scratch/y95/kgagalova/software/EVidenceModeler-v2.1.0" \

# make sure that the genomes and rnareads both have the same name denoted by the wildcard. In this case, the reads should have their prefix before _1.fastq.gz and _2.fastq.gz which should match with the name of the genomes.
#--funannotate_db "/scratch/y95/kgagalova/test_mycomanagement/mycomanagement_pilot/work/funannotate_db" \
