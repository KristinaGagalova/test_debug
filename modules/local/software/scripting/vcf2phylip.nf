process VCF_TO_PHYLIP {

    tag "${meta.id}"
    label 'small_task'
    label 'python'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/python:3.12':
        null }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.fasta"), optional: true, emit: fasta
    tuple val(meta), path("*.phy"), optional: true, emit: phylip
    tuple val(meta), path("*.nex"), optional: true, emit: nexus
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_name = params.refname ?: meta.id
    def output_fasta = task.ext.output_fasta ? "--fasta" : ""
    def output_phylip = task.ext.output_phylip ? "--phylip-sequential" : ""
    def output_nexus = task.ext.output_nexus ? "--nexus" : ""
    def prefix = task.ext.prefix ?: ""
    """
    python3 vcf2phylip.py \\
        -i ${vcf} \\
        ${output_fasta} \\
        ${output_phylip} \\
        ${output_nexus} \\
        ${prefix ? "--output-prefix ${prefix}" : ""} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        vcf2phylip: \$(python3 vcf2phylip.py --version 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def ref_name = params.refname ?: meta.id
    def output_fasta = task.ext.output_fasta ? "--fasta" : ""
    def output_nexus = task.ext.output_nexus ? "--nexus" : ""
    def prefix = task.ext.prefix ?: ref_name
    """
    # Create output files based on requested formats
    if [ "${output_fasta}" != "" ]; then
        touch ${prefix}.fasta
    fi
    if [ "${output_phylip}" != "" ]; then
        touch ${prefix}.phy
    fi
    if [ "${output_nexus}" != "" ]; then
        touch ${prefix}.nex
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2 || echo "3.12.0")
        vcf2phylip: "unknown"
    END_VERSIONS
    """
}