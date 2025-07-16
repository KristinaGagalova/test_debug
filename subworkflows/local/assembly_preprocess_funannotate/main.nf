//
// Clean and filter assembly
//

include { GAAS_FASTACLEANER }                                   from '../../../modules/local/software/gaas/fastacleaner'
include { GAAS_FASTASTATISTICS }                                from '../../../modules/local/software/gaas/fastastatistics'
include { GAAS_FASTAFILTERBYSIZE as GAAS_ASSEMBLYFILTERBYSIZE } from '../../../modules/local/software/gaas/fastafilterbysize'
include { FUNANNOTATE_CLEAN as FUNANNOTATE_CLEAN_GEN }          from '../../../modules/local/software/funannotate/funannotateclean.nf'
include { FUNANNOTATE_SORT }                                    from '../../../modules/local/software/funannotate/funannotatesort.nf'


workflow ASSEMBLY_PREPROCESS_FUNANNOTATE {

    take:
    genome // tuple: [meta, path]

    main:
    
    fasta_filt_clean        = Channel.empty()
    fasta_filt_clean_sort   = Channel.empty()
    fasta_stats             = Channel.empty()

    FUNANNOTATE_CLEAN_GEN(genome)
		.fasta
		.set { fasta_filt_clean }
    
    FUNANNOTATE_SORT(fasta_filt_clean)
                .fasta
                .set { fasta_filt_clean_sort }

    GAAS_FASTASTATISTICS(fasta_filt_clean)
		.stats
		.set { fasta_stats }

    emit:
    fasta_filt_clean_sort
    fasta_stats
}
