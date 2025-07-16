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
#export PATH=/scratch/y95/pmisiun/signalp-4.1:$PATH
#export PATH=/scratch/y95/kgagalova/src/signalp-5.0b/bin:$PATH

nextflow run ./main.nf -profile pawsey_setonix,singularity \
		--tool annotate-pangenome \
		--genomes "/scratch/y95/kgagalova/mp-cleanup/genomes_soft/*.fa" \
		--transcripts "/scratch/y95/kgagalova/mp-cleanup/transcripts/*fasta" \
		--gmark_db "/scratch/y95/kgagalova/src/genemark/gmes_linux_64_4/" \
		--funannotate_db "/scratch/y95/kgagalova/mp-cleanup/mycoprocessor-final/work/funannotate_db" \
		--cqpm true \
		--signalp_path "/scratch/y95/kgagalova/mp-cleanup/test/signalp-4.1g-work/signalp-4.1g-install/" \
		--enable_conda false \
		-with-dag flowchart.png \
		-resume
		# --augustus_config "/scratch/y95/kgagalova/software/Augustus/config" \
        # --evm_home "/scratch/y95/kgagalova/software/EVidenceModeler-v2.1.0" \
