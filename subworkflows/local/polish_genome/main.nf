//
// Subworkflow - polish genome
//

include { BWA_INDEX }                            from '../../../modules/local/software/bwa/indexgenome'
include { BWAMEM_MAPALL as BWAMEM_R1 }           from '../../../modules/local/software/bwa/mapallreads'
include { BWAMEM_MAPALL as BWAMEM_R2 }           from '../../../modules/local/software/bwa/mapallreads'
include { POLYPOLISHFILT }                       from '../../../modules/local/software/polypolish/polishgenome'
include { POLYPOLISHPOLISH }                     from '../../../modules/local/software/polypolish/polishgenome'

workflow POLISH_GENOME {

    take:
    sample        // tuple [ sample, reads_pair1, reads_pair2, long_reads, contamination_bool, assembly ]


    main:
    sample_pair_genome = sample
        .map { meta, reads1, reads2, long_reads, contam_bool, genome ->
        [meta, reads1, reads2]
        }

    genome = sample
        .map { meta, reads1, reads2, long_reads, contam_bool, genome ->
        [meta, genome]
        }
    

	// Pass the genome to the BWA_INDEX process
	indexed_genome = BWA_INDEX(genome)
	
        // Create channels for each R1 and R2 plus genome
	sample_pair_genome.join(indexed_genome.index)
		.map { sample, reads1, reads2, genome -> tuple(sample, reads1, genome)}
		.set { reads_gen_r1 }

	sample_pair_genome.join(indexed_genome.index)
                .map { sample, reads1, reads2, genome -> tuple(sample, reads2, genome)}
                .set { reads_gen_r2 }	

	// Map reads individually using the indexed genome output as input for the index
        reads_R1 = BWAMEM_R1(
		reads_gen_r1,
        	"R1"
		)

         reads_R2 = BWAMEM_R2(
		reads_gen_r2,
                "R2"
            )

        // Pair reads together after mapping
        reads_R1.reads_mapped.join(reads_R2.reads_mapped)
                .set { ch_readsR1_readsR2 }

        // Polypolish filtering step
        FILT_READS = POLYPOLISHFILT(ch_readsR1_readsR2)

	// Pair filtered reads and genome
	FILT_READS.reads_filt.join(genome)
		.set { ch_readsR1_readsR2_genome }

        // Polypolish polishing genome step
        polished_genome = POLYPOLISHPOLISH(
                ch_readsR1_readsR2_genome
        	)

    emit:
	filt_reads  = FILT_READS.reads_filt
    assembly = polished_genome.genome_polished
}