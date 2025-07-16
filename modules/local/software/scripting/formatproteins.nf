process FORMAT_PROTEINS {
    label 'small_task'

    input:
    path(proteins_cds)

    output:
    tuple path("00_proteomes"), path("00_transcripts"), emit: formatted_proteins_cds

    script:
    """
    mkdir 00_proteomes
    mkdir 00_transcripts

    mv *.proteins.fa 00_proteomes
    mv *.cds.fa 00_transcripts

    for file in 00_proteomes/*; do
        mv "\$file" "\${file/.proteins/}"
        temp_file="\${file}.tmp"
        awk '/^>/ {if (seqlen) {print seq}; print; seq=""; seqlen=0; next} {seq=seq""\$0; seqlen+=length(\$0)} END {if (seqlen) print seq}' "\${file/.proteins/}" > "\$temp_file"
        mv "\$temp_file" "\${file/.proteins/}"
        sed -i -E 's/PGN\\.([0-9]+)/CQPM_\\1 CQPM_\\1/' "\${file/.proteins/}"
    done

    for file in 00_transcripts/*; do
        mv "\$file" "\${file/.cds/}"
        temp_file="\${file}.tmp"
        awk '/^>/ {if (seqlen) {print seq}; print; seq=""; seqlen=0; next} {seq=seq""\$0; seqlen+=length(\$0)} END {if (seqlen) print seq}' "\${file/.cds/}" > "\$temp_file"
        mv "\$temp_file" "\${file/.cds/}"
        sed -i -E 's/PGN\\.([0-9]+)/CQPM_\\1 CQPM_\\1/' "\${file/.cds/}"
    done
    """
}