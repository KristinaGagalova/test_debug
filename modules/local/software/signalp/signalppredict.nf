process SIGNALP_RUN {

    label 'medium_task'
    tag "${meta}"

    // software must be installed locally
    input:
    tuple val(meta), path(chunk_file)

    output:
    tuple val(meta), path("CQ_${meta}_${chunk_file}"), emit: signalp_output

    script:
    """
    ${params.signalp_path}/signalp ${chunk_file} > CQ_${meta}_${chunk_file}
    """
}

process SIGNALP_AGGREGATE {

    input:
    tuple val(meta), path(signalp_outputs) // signal outputs contains file chunks

    output:
    tuple val(meta), path("Secretome_${meta}.txt"), emit: secretome

    script:
    """
    cat ${signalp_outputs} | grep -v "#" | awk '(\$10 == "Y"){print \$1" "\$5}' > Secretome_${meta}.txt
    """
}
