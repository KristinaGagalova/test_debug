//
// Subworkflow structures and ideas from https://github.com/nf-core/genomeannotator
// https://github.com/ERGA-consortium/pipelines/blob/main/annotation/nextflow/modules/local/funannotate/FunannotatePredict.nf
//

include { GAAS_FASTACLEANER }                               from '../../../modules/local/software/gaas/fastacleaner'
include { EXONERATE_FASTACLEAN }                            from '../../../modules/local/software/exonerate/fastaclean'
include { SEQKIT_DEDUP_FASTA }                              from '../../../modules/local/software/seqkit/removeduplicated'
include { FUNANNOTATEPASA_GENE_PREDICTION_BAM }             from '../../../modules/local/software/funannotate/funannotatepredict'
include { FUNANNOTATEPASA_GENE_PREDICTION_CQ }              from '../../../modules/local/software/funannotate/funannotatepredict'
include { PASA_SEQCLEAN }                                   from '../../../modules/local/software/pasa/seqclean'
include { PASA_ALIGNASSEMBLE }                              from '../../../modules/local/software/pasa/alignassemble'
include { PASA_ASMBLSTOTRAINING }                           from '../../../modules/local/software/pasa/assemblstotraining'
include { HELPER_PASA2TRAINING }                            from '../../../modules/local/software/pasa/pasa2training'
include { CQ_TRANSCRIPTOME_PREDICT }                        from '../../../modules/local/software/codingquarry/cqtranscriptomepredict'
include { CQ_PATHOGEN_MODE }                                from '../../../subworkflows/local/run_codingquarry-pm/main'
include { FUNANNOTATE_GENE_PREDICTION_TRAINED }             from '../../../modules/local/software/funannotate/funannotatepredict'
include { CONCAT_ANNOTATIONS as CONCAT_ANNOTATION_PROTSEQ } from '../../../modules/local/software/scripting/concatenatefasta'
include { CONCAT_ANNOTATIONS as CONCAT_ANNOTATION_CDS }     from '../../../modules/local/software/scripting/concatenatefasta'


workflow FUNANNOTATEPREDICT_REF_PIPELINE {

    take:
    genome_training // tuple val(meta), path(fasta), path(training), path(bam)
    proteins // file path (may be null)
    funannotate_setup

    main:    
    genome_training
		.map { vector ->
    		// Drop the last element of the vector
    		vector[0..-3]}
		.set { genome }


    // run gene prediction module, using CQ approach
    annotated = FUNANNOTATE_GENE_PREDICTION_TRAINED(
        funannotate_setup,
        genome_training,
        proteins,
        params.species
        )
    
    protseq = annotated.protseq
    cds = annotated.cds


    genome_annotations = genome.join(annotated.gff)

    if (params.cqpm) {

            annotated_cqpm = CQ_PATHOGEN_MODE(genome_annotations).output
            unmerged_protseq = annotated.protseq.join(CQ_PATHOGEN_MODE.out.protseq)
            unmerged_cds = annotated.cds.join(CQ_PATHOGEN_MODE.out.cds)
            protseq = CONCAT_ANNOTATION_PROTSEQ(unmerged_protseq)
            cds = CONCAT_ANNOTATION_CDS(unmerged_cds)

    } else {

       annotated_cqpm = null

    }

    
    emit:
       funannotate_predictions = annotated.output
       annotated_cqpm = annotated_cqpm
       protein_fasta = protseq
       cds_fasta = cds
       //versions = annotated.versions
       
}



workflow FUNANNOTATEPREDICT_ALT_PIPELINE {

    take:
    genome // tuple val(meta), path(fasta)
    genome_bam
    proteins // file path (may be null)
    transcripts // file path (may be null)
    funannotate_setup

    main:    
    // process transcripts if any
    if (transcripts) {

        GAAS_FASTACLEANER(
            transcripts
        )

        EXONERATE_FASTACLEAN(
            GAAS_FASTACLEANER.out.fasta
        )

        SEQKIT_DEDUP_FASTA(
            EXONERATE_FASTACLEAN.out.fasta
        )

        // run pasa predictions
        PASA_SEQCLEAN(
            SEQKIT_DEDUP_FASTA.out.fasta
        )

        PASA_ALIGNASSEMBLE(
            genome,
            PASA_SEQCLEAN.out.fasta,
            params.config_pasa,
            params.config_pasa_dir,
            params.max_intron_size
        )
 
        PASA_ASMBLSTOTRAINING(
            PASA_ALIGNASSEMBLE.out.pasa_out
        )


	ch_genome_pasa = genome.join(PASA_ASMBLSTOTRAINING.out.gff)
        CQ_TRANSCRIPTOME_PREDICT(
        ch_genome_pasa
        )

        	// run gene prediction module, using CQ approach
	annotated = FUNANNOTATEPASA_GENE_PREDICTION_CQ(
                funannotate_setup,
                CQ_TRANSCRIPTOME_PREDICT.out, //This has been edited to be a larger tuple including CQ gff
               	proteins,
                transcripts,
                params.species
                )

        protseq = annotated.protseq
        cds = annotated.cds

	}

    genome_annotations = genome.join(annotated.gff)
    if (params.cqpm) {

            annotated_cqpm = CQ_PATHOGEN_MODE(genome_annotations).output
            unmerged_protseq = annotated.protseq.join(CQ_PATHOGEN_MODE.out.protseq)
            unmerged_cds = annotated.cds.join(CQ_PATHOGEN_MODE.out.cds)
            protseq = CONCAT_ANNOTATION_PROTSEQ(unmerged_protseq)
            cds = CONCAT_ANNOTATION_CDS(unmerged_cds)

    } else {

       annotated_cqpm = null

    }

    
    emit:
       pasa_predictions = PASA_ASMBLSTOTRAINING.out.gff
       funannotate_predictions = annotated.output
       annotated_cqpm = annotated_cqpm
       protein_fasta = protseq
       cds_fasta = cds
       //versions = annotated.versions
       
}
