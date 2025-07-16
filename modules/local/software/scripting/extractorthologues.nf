process EXTRACT_ORTHOLOGUES {
    label 'medium_task'

    input:
    path(protein_ortho)
    tuple path(prots), path(cds)
    val(prefix)

    output:
    path("03_parsed_orthogroups/representatives.prot.fasta"),   emit: protein_reps
    path("03_parsed_orthogroups/representatives.cds.fasta"),    emit: cds_reps

    script:
    """
    mkdir 02_raw_orthogroups
    mkdir 03_parsed_orthogroups

    ##extract raw ortholog groups
    cat ${protein_ortho} | extract_orthogroups_fasta.pl ${prots} 02_raw_orthogroups/ ${prefix}

    ##Select 1 representative sequence for each group
    select_orthogroups_reps.pl 02_raw_orthogroups ${cds} 03_parsed_orthogroups ${prefix}
    """
}