include { SHORT_READ_CONTAM_FILTERING }         from '../../local/contamination_filtering/main'

include { TRIM_READS }                          from '../../../modules/local/software/fastp/trimreads'
include { BBMERGE }                             from '../../../modules/local/software/bbmap/bbmerge'
include { SPADES_ASSEMBLY }                     from '../../../modules/local/software/spades/spadesassembly'
include { MITOZ_ASSEMBLY }                      from '../../../modules/local/software/mitoz/mitozassembly'
include { ASSEMBLY_STATS}                       from '../../../modules/local/software/scripting/assemblystats'
include { CONCAT_FASTQ }                        from '../../../modules/local/software/scripting/concatenatefastq'
include { BLAST_MITOGENOME }                    from '../../../modules/local/software/blast/blastmitogenome'
include { SAMTOOLS_INDEX_GENOME }               from '../../../modules/local/software/samtools/indexgenome'
include { EXCLUDE_MITOGENOME }                  from '../../../modules/local/software/scripting/excludemitogenome'
include { FINAL_SHORT_READ_CLEANUP}             from '../../../modules/local/software/scripting/finalshortreadcleanup'
include { FINAL_SHORT_READ_CLEANUP_NO_MITOZ}    from '../../../modules/local/software/scripting/finalshortreadcleanup'

workflow ASSEMBLE_SHORT_READS {
    take:
    reads

    main:
    //remove any null fields (should only be long reads)
    reads.map { tuple -> tuple.findAll { it != null } }.set { reads }

    //create channel for just reads
    reads.map { name, reads1, reads2, contam -> tuple(name, reads1, reads2) }.set { only_reads }

    //create channel for contamination filtering T/F
    reads.map { name, reads1, reads2, contam -> tuple(name, contam) }.set { contam_bool }


    only_reads.branch {
    multi: it[1] instanceof Collection  // Branch if the second element is a List
    single: !(it[1] instanceof Collection)  // Branch if the second element is not a List
    other: true
    }
    .set{ branched_reads }


    CONCAT_FASTQ(branched_reads.multi)
    .set { concatenated_multis }

    concatenated_multis.mix(branched_reads.single)
    .set { mixed_paired_reads }

    TRIM_READS(mixed_paired_reads)
    .set { trimmed_reads }

    //add contamination T/F back
    trimmed_reads.join(contam_bool)
    .set { trimmed_with_contam }

    trimmed_with_contam.view()


    trimmed_with_contam.branch {
    contam_true: it[3] == "True" // Contam = True
    contam_false: it[3] == "False"  // Contam = False
    }
    .set { contamination_channels }


    //send samples with contam = True through contam filtering pipeline
    SHORT_READ_CONTAM_FILTERING(contamination_channels.contam_true)
    .set{ contamination_filtered }

    //remove contam feild from contam=F samples
    contamination_channels.contam_false.map { name, reads1, reads2, contam -> tuple(name, reads1, reads2) }.set { contam_unfiltered_reads }

    //join contamination filtered and unfiltered channels
    contamination_filtered.contam_filtered_reads.mix(contam_unfiltered_reads).
    set { final_reads_ch }


    //read merging
    BBMERGE(final_reads_ch)
    .set { merged_reads }

    //spades assembly
    SPADES_ASSEMBLY(merged_reads.reads)
    .set {assemblies}

    //index genome
    SAMTOOLS_INDEX_GENOME(assemblies.assembly)
    .set { indices }

    //mitoz assembly
    MITOZ_ASSEMBLY(trimmed_reads)
    .set{ mitoz_out }

    assemblies.assembly.join(mitoz_out, remainder: true)
    .set { mito_and_ncl }


    mito_and_ncl.branch {
    mitoz_pass:  it[2] != null 
    mitoz_fail: it[2] == null 
    }
    .set { mito_and_ncl_branched }

    mito_and_ncl_branched.mitoz_fail.map { tuple -> tuple.findAll { it != null } }.set { mito_and_ncl_mitoz_fail }
    
    mito_and_ncl_mitoz_fail.view()

    FINAL_SHORT_READ_CLEANUP_NO_MITOZ(mito_and_ncl_mitoz_fail)
    .set{ final_assemblies_mitoz_fail }

    //generate assembly stats
    // ASSEMBLY_STATS(mito_and_ncl.mitoz_pass)
    // .set { assembly_stats }


    // //Cleanup steps
    BLAST_MITOGENOME(mito_and_ncl_branched.mitoz_pass)
    .set { mito_blast_results }


    mito_blast_results.blast_table.join(assemblies.assembly)
    .set { blast_genome }

    blast_genome.join(indices.index)
    .set { blast_genome_index }

    EXCLUDE_MITOGENOME(blast_genome_index)
    .set { mitogenome_excluded }

    mitogenome_excluded.excluded_fasta.join(mitoz_out)
    .set { ncl_and_mito }

    FINAL_SHORT_READ_CLEANUP(ncl_and_mito)
    .set { final_assemblies_mitoz_pass }

    final_assemblies_mitoz_pass.final_assembly.mix(final_assemblies_mitoz_fail.final_assembly)
    .set { final_assemblies_all }

    final_assemblies_all.view()

    
}