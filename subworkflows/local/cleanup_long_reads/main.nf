include { SEQKIT_LISTRENAME } from '../../../modules/local/software/seqkit/renamefasta'
include { SEQKIT_RENAMEGEN } from '../../../modules/local/software/seqkit/renamefasta'
include { MITO_ALIGN } from '../../../modules/local/software/minimap/check_mito'
include { SEQKIT_MITOTAG } from '../../../modules/local/software/seqkit/renamefasta'

workflow CLEANUP_GENOME {

    take:
    genome        // tuple [ sample, genome ]
	mito_ref      // path to download

    main:


	LIST_NAMS = SEQKIT_LISTRENAME(genome)

    LIST_NAMS.list_scafs.join(genome)
		.set { ch_genome_orig_nams }

    ASSEMBLY_RENAM = SEQKIT_RENAMEGEN(ch_genome_orig_nams)

    LIST_MITO = MITO_ALIGN(ASSEMBLY_RENAM.renamed_fasta, mito_ref)

    LIST_MITO.list_mito.join(ASSEMBLY_RENAM.renamed_fasta)
		.set { ch_renamed_fasta_mito }

    FINAL_GENOME = SEQKIT_MITOTAG(ch_renamed_fasta_mito)

    emit:
    assembly = FINAL_GENOME.mito_fasta

}