include { FUNANNOTATEPREDICT_ALT_PIPELINE }  from '../subworkflows/local/run_annotate-pangenome/main'
include { FUNANNOTATEPREDICT_REF_PIPELINE }  from '../subworkflows/local/run_annotate-pangenome/main'
include { ASSEMBLY_PREPROCESS_FUNANNOTATE }  from '../subworkflows/local/assembly_preprocess_funannotate/main'

include { FUNANNOTATE_SETUP }                from '../modules/local/software/funannotate/funannotatesetup'
include { CONCAT_FASTA }                     from '../modules/local/software/scripting/concatenatefasta'
include { FUNANNOTATE_TRAIN }                from '../modules/local/software/funannotate/funannotatetrain'
include { PROTEIN_ORTHO_CLUSTER }            from '../modules/local/software/proteinortho/proteinorthoclustering'
include { PREDICT_REPRESENTATIVES }			 from '../modules/local/software/metaeuk/predictrepresentatives'
include { REVISE_ORTHOGROUPS }               from '../modules/local/software/scripting/reviseorthogroups'
include { FORMAT_READ_HEADERS }				 from '../modules/local/software/scripting/formatreadheaders'

if (params.genomes) { ch_genomes = Channel.fromPath(params.genomes, checkIfExists: true) } else { exit 1, 'No assembly genomes specified!' }
if (params.transcripts) { transcripts = Channel.fromPath(params.transcripts).collect() } else { transcripts = "EMPTY" }
if (params.proteins) {proteins = Channel.fromPath(params.proteins) } else { proteins = "EMPTY" }
if (params.rnareads) {ch_rnareads = Channel.fromFilePairs(params.rnareads, checkIfExists: true)} else { ch_rnareads = "EMPTY" }


if (params.gmark_db == null ) exit 1, "Parameter 'gmark_db' is null. Exiting the run. Please provide path the genemark gmes_petap.pl"
//if (params.funannotate_db == null ) exit 1, "Parameter 'funannotate_db' is null. Exiting the run. Please provide a valid path to the DB"

if (params.cqpm && params.signalp_path == null) exit 1, "Please download and register signalP v4.1: https://services.healthtech.dtu.dk/services/SignalP-4.1/"

//------------------------ FUNCTIONS

def processChannels(ch_input) {
    return ch_input.map { path ->
        def name = "${path.baseName}"
        tuple(name, path)
    }
}

workflow ANNOTATE_PANGENOME_WORKFLOW {

	//Establish genome channel and preprocess
	def ch_genomesNam_unprocessed = processChannels(ch_genomes)
        ch_genomesNam = ASSEMBLY_PREPROCESS_FUNANNOTATE(ch_genomesNam_unprocessed).fasta_filt_clean_sort

	FUNANNOTATE_SETUP(
    	    params.buscodb,
    	    "${baseDir}/work/funannotate_db"	
    	    )
        params.funannotate_db = "${baseDir}/work/funannotate_db"

	if (ch_rnareads != "EMPTY") {

		// FORMAT_READ_HEADERS(ch_rnareads)
		// .set {ch_rnareads_formatted}
		
		// filter out non matching elements
		ch_genomesNam
			.join(ch_rnareads, remainder: true)
			.set { ch_genomes_rnareads }

		ch_genomes_rnareads
    			.ifEmpty {exit 1, 
       			"There are no common files between the reads and the genomes. Please check your input. Stopping execution."
    		}

		//split genome channel into those with and without rna
		ch_genomes_rnareads.branch {
    		ch_genomes_with_rna: it[-1] != null
    		ch_genomes_without_rna: it[-1] == null
		}
		.set{ ch_genomes_rnareads_split }
		
		//define split channels
		ch_genomes_with_rna = ch_genomes_rnareads_split.ch_genomes_with_rna
		ch_genomes_rnareads_split.ch_genomes_without_rna
		.map { vector ->
    	 	// Drop the last element of the vector
    	 	vector[0..-2]}
		.set{ ch_genomes_without_rna }
		
		//run funannotate train for those with rnaseq
		FUNANNOTATE_TRAIN(
			ch_genomes_with_rna,
			params.species
		)

        ch_genomes_with_rna
		.map { vector ->
    		// Drop the last element of the vector
    		vector[0..-2]}
		.set { ch_genomes_with_rna_trained }
		

		// create tuples for input	
		ch_genome_training = ch_genomes_with_rna_trained
			.join(FUNANNOTATE_TRAIN.out.training)
		
		ch_genome_training_bam = ch_genome_training
			.join(FUNANNOTATE_TRAIN.out.bam)
        
		// run genomes with rna training through ref pipeline
		FUNANNOTATEPREDICT_REF_PIPELINE(
			ch_genome_training_bam,
			proteins,
			FUNANNOTATE_SETUP.out.versions
		)

		//run those without rna training through alt pipeline
		FUNANNOTATEPREDICT_ALT_PIPELINE(
		ch_genomes_without_rna,
		"EMPTY",
		proteins,
		CONCAT_FASTA(transcripts.mix(FUNANNOTATE_TRAIN.out.transcripts).collect()).merged, //This will also include any transcripts derived from ref genomes with funannotate train
		FUNANNOTATE_SETUP.out.versions
	 	)

		protseq_meta = FUNANNOTATEPREDICT_REF_PIPELINE.out.protein_fasta.mix(FUNANNOTATEPREDICT_ALT_PIPELINE.out.protein_fasta)
		cds_meta = FUNANNOTATEPREDICT_REF_PIPELINE.out.cds_fasta.mix(FUNANNOTATEPREDICT_ALT_PIPELINE.out.cds_fasta)
		cds_protseq_meta = protseq_meta.join(cds_meta)

		cds_protseq_meta
		.map { vector ->
    		// Drop the first element of the vector
    		vector[1..-1]}
		.collect()
		.set { cds_protseq }
		



	
	} else {

		//When no genomes have rna reads, run them all through alt pipeline
		ch_genomes_bam = ch_rnareads
		FUNANNOTATEPREDICT_ALT_PIPELINE(
		ASSEMBLY_PREPROCESS_FUNANNOTATE.out.fasta_filt_clean_sort,
		ch_genomes_bam,
		proteins,
		CONCAT_FASTA(transcripts).merged,
		FUNANNOTATE_SETUP.out.versions
	 	)

		

		protseq_meta = FUNANNOTATEPREDICT_ALT_PIPELINE.out.protein_fasta
		cds_meta = FUNANNOTATEPREDICT_ALT_PIPELINE.out.cds_fasta
		cds_protseq_meta = protseq_meta.join(cds_meta)

		cds_protseq_meta
		.map { vector ->
    		// Drop the first element of the vector
    		vector[1..-1]}
		.collect()
		.set { cds_protseq }

		
	}



}
