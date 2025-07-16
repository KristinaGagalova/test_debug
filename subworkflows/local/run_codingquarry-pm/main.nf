include { CQ_EXTRACT_CDS }           from '../../../modules/local/software/codingquarry/cqtranscriptomepredict'
include { CQ_CDS_FILTER_PROCESS }    from '../../../modules/local/software/codingquarry/cqtranscriptomepredict'
include { SIGNALP_RUN }              from '../../../modules/local/software/signalp/signalppredict'
include { SIGNALP_AGGREGATE }        from '../../../modules/local/software/signalp/signalppredict'
include { CQ_PATHOGENMODE }          from '../../../modules/local/software/codingquarry/cqtranscriptomepredict'            

//----------------------------------------------

workflow CQ_PATHOGEN_MODE {

    take:    
    input_fasta // tuple name, path(fasta_genome), path(gff)

    main:
    CQ_EXTRACT_CDS(input_fasta)
    CQ_PROCESSED = CQ_CDS_FILTER_PROCESS(CQ_EXTRACT_CDS.out.predicted_cds)
   
    CQ_PROCESSED.split_files
    	.flatMap { name, filePaths -> filePaths.collect { filePath -> [name, filePath] } }
    	.set { files_to_process }
     
    SIGNALP_RUN(files_to_process).signalp_output
    		.groupTuple()
    		.set { grouped_files_output }
    // collect all the output
    SIGNALP_AGGREGATE(grouped_files_output)
    SIGNALP_AGGREGATE.out.secretome.view()
    joined_secretome = input_fasta.join(SIGNALP_AGGREGATE.out.secretome)    
    CQ_PATHOGENMODE(joined_secretome)

    emit:
    output = CQ_PATHOGENMODE.out.output
    cds = CQ_PATHOGENMODE.out.cds
    protseq = CQ_PATHOGENMODE.out.protseq
}
