#!/usr/bin/env nextflow

//Call DSL2
nextflow.enable.dsl=2

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
//ver = manifest.version

//--------------------------------------------------------------------------------------------------------
// Validation - validation from nf-core for input results

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

if (params.help) {
   log.info paramsHelp("nextflow run ...")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

//--------------------------------------------------------------------------------------------------------

// This part calls the workflows
workflow_input = params.tool
switch (workflow_input) {
    case ["occultercut"]:
        include { OCCULTERCUT_WORKFLOW } from './workflows/run_occultercut.nf'
	break;
    case ["funannotate"]:
	include { FUNANNOTATE_WORKFLOW } from './workflows/run_funannotate.nf'
	break;
	case ["annotate-pangenome"]:
	include { ANNOTATE_PANGENOME_WORKFLOW } from './workflows/run_annotate-pangenome.nf'
	break;
	case ["assemble-pangenome"]:
	include { ASSEMBLE_PANGENOME_WORKFLOW } from './workflows/run_assemble_pangenome.nf'
	break;
	case ["orthology-clustering"]:
	include { ORTHOLOGY_CLUSTERING_WORKFLOW } from './workflows/run_orthology_clustering.nf'
	break;
	case ["snp-phylo"]:
	include { SNPS_PHYLO_WORKFLOW } from './workflows/run_snps_phylo.nf'
	break;

}


// Main workflow used to select from themes and tools
workflow {

	/*
	* WORKFLOW THEME #: add name and description here 
	*/
	if (params.tool == "occultercut") {
		OCCULTERCUT_WORKFLOW()

	} else if (params.tool == "annotate-pangenome") {
		ANNOTATE_PANGENOME_WORKFLOW()
	
	} else if (params.tool == "assemble-pangenome") {
		ASSEMBLE_PANGENOME_WORKFLOW()

	} else if (params.tool == "orthology-clustering") {
		ORTHOLOGY_CLUSTERING_WORKFLOW()

	} else if (params.tool == "snp-phylo") {
		SNPS_PHYLO_WORKFLOW()

	} else {
	
		println("Please provide the correct input options")

	}		 
}
