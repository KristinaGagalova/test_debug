include { POLISH_GENOME } from '../../local/polish_genome/main'
include { CLEANUP_GENOME } from '../../local/cleanup_long_reads/main'

include { FLYE_ASSEMBLE_MOD } from '../../../modules/local/software/flye/flyeassembly'
include { NANO_COMPARE_FASTA } from '../../../modules/local/software/nanocompare/nanocomparefasta'
include { NANO_COMPARE_FASTQ } from '../../../modules/local/software/nanocompare/nanocomparefastq'
include { MITOCONDRION_DOWNLOAD } from '../../../modules/local/software/download/downloadmito'


workflow ASSEMBLE_NANO {

    take:
    sample // tuple[meta, reads1, reads2, long_reads, contam_bool]

    main:
    reads = sample
        .map { meta, reads1, reads2, long_reads, contam_bool ->
        [meta, long_reads]
        }
    

    only_fastq = reads
        .map {meta, reads ->
        [reads[0]]
        }
        .collect()

    only_fasta = reads
        .map {meta, reads ->
        [reads[1]]
        }
        .collect()

    //QC
    NANO_COMPARE_FASTQ(only_fastq)
    NANO_COMPARE_FASTA(only_fasta)
    
    //ASSEMBLY
    FLYE_ASSEMBLE_MOD(reads)

    sample.join(FLYE_ASSEMBLE_MOD.out.fasta)
    .set { samples_with_genome }

    samples_with_genome.branch {
        only_long_reads: it[1] == null && it[2] == null && it[3] != null
        paired_and_long_reads: it[1] != null && it[2] != null && it[3] != null
    }
    .set {branched_assembled_long_reads}

    //POLISHING FOR SAMPLES WITH SHORT READS
    polished_assemblies = POLISH_GENOME(branched_assembled_long_reads.paired_and_long_reads)

    only_long_reads = branched_assembled_long_reads.only_long_reads.map { 
        meta, reads1, reads2, long_reads, contam_bool, genome ->
        [meta, genome]
        }

    all_long_read_assemblies = only_long_reads.mix(polished_assemblies.assembly)

    all_long_read_assemblies.view()

    //DOWNLOAD MITO
    MITO_CHECK = MITOCONDRION_DOWNLOAD(params.mito_db)

    //FINAL CLEANUP
    CLEANED_GENOME = CLEANUP_GENOME(all_long_read_assemblies, MITO_CHECK.mito_ref)

    emit:
    assembly = CLEANED_GENOME.assembly
}