include { FORMAT_PROTEINS }                  from '../modules/local/software/scripting/formatproteins'
include { PROTEIN_ORTHO_CLUSTER }            from '../modules/local/software/proteinortho/proteinorthoclustering'
include { EXTRACT_ORTHOLOGUES }              from '../modules/local/software/scripting/extractorthologues'
include { PREDICT_REPRESENTATIVES }			 from '../modules/local/software/metaeuk/predictrepresentatives'
include { GENE_TABLES }                      from '../modules/local/software/scripting/reviseorthogroups'
include { REVISE_ORTHOGROUPS }               from '../modules/local/software/scripting/reviseorthogroups'
include { REVISE_ORTHOGROUPS_ALT }           from '../modules/local/software/scripting/reviseorthogroups'

if (params.genomes) { ch_genomes = Channel.fromPath(params.genomes, checkIfExists: true) } else { exit 1, 'No assembly genomes specified!' }
if (params.cds) { ch_cds = Channel.fromPath(params.cds, checkIfExists: true) } else { println 'No cds specified!' }
if (params.prot) { ch_prot = Channel.fromPath(params.prot, checkIfExists: true) } else { println 'No proteins specified!' }
def pangene_inputs_given = params.pangene_cds && params.pangene_prot
if (pangene_inputs_given) {
    ch_pangene_cds = file(params.pangene_cds)
    ch_pangene_prot = file(params.pangene_prot)
} else {
    log.warn 'Missing pangene input files! Workflow will not proceed down this path.'
}


def processChannels(ch_input, suffixToRemove) {
    return ch_input.map { path ->
        def name = "${path.baseName}".replace(suffixToRemove, "")
        tuple(name, path)
    }
}

workflow ORTHOLOGY_CLUSTERING_WORKFLOW {
    def ch_genomes_processed = processChannels(ch_genomes, "")

    if (pangene_inputs_given) {
        println "this path"
        pangene_prot = ch_pangene_prot
        pangene_cds = ch_pangene_cds

        PREDICT_REPRESENTATIVES(ch_genomes_processed,
                            pangene_prot)

        GENE_TABLES(PREDICT_REPRESENTATIVES.out.metaeuk_output_all,
                pangene_prot,
				params.prefix)

        REVISE_ORTHOGROUPS_ALT(PREDICT_REPRESENTATIVES.out.metaeuk_output_all.collect(), 
                            GENE_TABLES.out.pav_table.collect(), 
                            pangene_prot,
                            pangene_cds,
							params.prefix,
							params.cull_groups)
    } 
    else {
        println "other path"
        def ch_cds_processed = processChannels(ch_cds, ".cds")
        def ch_prot_processed = processChannels(ch_prot, ".proteins")

        cds_protseq_meta = ch_cds_processed.join(ch_prot_processed)

        cds_protseq_meta
            .map { vector ->
                // Drop the first element of the vector
                vector[1..-1]}
            .collect()
            .set { cds_protseq }

        cds_protseq.view()
        
        //format proteins for protein ortho
        FORMAT_PROTEINS(cds_protseq)

        //protein ortho clustering
        PROTEIN_ORTHO_CLUSTER(FORMAT_PROTEINS.out.formatted_proteins_cds)

        //extract orthologues and representatives
        EXTRACT_ORTHOLOGUES(PROTEIN_ORTHO_CLUSTER.out.protein_ortho,
                            FORMAT_PROTEINS.out.formatted_proteins_cds,
							params.prefix)
        
        pangene_prot = EXTRACT_ORTHOLOGUES.out.protein_reps
        pangene_cds = EXTRACT_ORTHOLOGUES.out.cds_reps

        //metaeuk scan
        PREDICT_REPRESENTATIVES(ch_genomes_processed,
                                pangene_prot)

        GENE_TABLES(PREDICT_REPRESENTATIVES.out.metaeuk_output_all,
                    pangene_prot,
					params.prefix)

        //cull empty orthologue groups
        REVISE_ORTHOGROUPS(PREDICT_REPRESENTATIVES.out.metaeuk_output_all.collect(), 
                            GENE_TABLES.out.pav_table.collect(), 
                            pangene_prot,
                            pangene_cds,
							params.prefix)
    }
    
}
