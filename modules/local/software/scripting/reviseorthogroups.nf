params.cull_groups = params.cull_groups ?: 0  // Set default to 0

process GENE_TABLES {
    label 'small_task'
    
    input:
    path(metaeuk_results)
    path(prot_reps)
    val(prefix)

    output:
    path("${metaeuk_results[0].baseName}.PAV_*.txt"), emit: pav_table

    script:
    """
    summary=${metaeuk_results[0].baseName}.PAV_counts.txt
    pav=${metaeuk_results[0].baseName}.PAV_table.txt

    #create list of orthogroups to search for
	list=\$(grep "^>"${prefix}"_" ${prot_reps} | sed 's/^>//')
    echo \$list

    for file in \$(find -name "*_metaeuk.fas" | sed 's/_metaeuk.fas//' | sed 's/\\.\\///'); do
        echo \$file
        headers=\$file"_metaeuk.headersMap.tsv"
        prot=\$file"_metaeuk.fas"
        nucl=\$file"_metaeuk.codon.fas"
        
        echo -n \$file >> \$summary
        echo -n \$file >> \$pav
        for gene in \$list; do    
            if [ -z "\$(grep "\$gene|" < \$headers)" ]; then
                matchlist="."
            else
                matchlist=\$(grep "\$gene|" < \$headers | awk '{print \$6}' | paste -sd ',')
            fi
            count=\$(grep -c "\$gene|" < \$headers) || true
            echo -n -e '\t' \$count >> \$summary
            echo -n -e '\t' \$matchlist >> \$pav
        done
        echo "" >> \$summary
        echo "" >> \$pav
    done
    """
}

process REVISE_ORTHOGROUPS { //generate representative ortholog group panel - this process is for initial run from scratch with supporting protein/transcriptome inputs
    label 'medium_task'

    input:
    path(metaeuk_results)
    path(pavs)
    path(prot_reps)
    path(cds_reps)
    val(prefix)

    output:
    path("revised_representatives.prot.fasta"), emit: prot_reps
    path("revised_representatives.cds.fasta"),  emit: cds_reps
    path("revised_PAV_counts.txt"), emit: orthogroup_PAV_counts
    path("revised_PAV_table.txt"), emit: orthogroup_PAV

    script:
    """
    reviseOrthogroups.sh ${prot_reps} ${cds_reps} ${prefix}
    """
}

process REVISE_ORTHOGROUPS_ALT { //map representative ortholog group panel - this process is for secondary runs where the panel has previously been generated and only needs to be mapped/annotated onto genome fasta input

    label 'medium_task'

    input:
    path(metaeuk_results)
    path(pavs)
    path(prot_reps)
    path(cds_reps)
    val(prefix)
	val(cull_groups)

    output:
    path("revised_representatives.prot.fasta"), emit: prot_reps
    path("revised_representatives.cds.fasta"),  emit: cds_reps
    path("revised_PAV_counts.txt"), emit: orthogroup_PAV_counts
    path("revised_PAV_table.txt"), emit: orthogroup_PAV

    script:
    """
    reviseOrthogroupsAlt.sh ${prot_reps} ${cds_reps} ${prefix} ${cull_groups}
    """
}