process GENERATE_TREE {

    tag "${meta.id}"
    label 'medium_task'
    label 'iqtree'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.4--h21ec9f0_0':
        null }"

    input:
    tuple val(meta), path(phylip)

    output:
    tuple val(meta), path("*.treefile")     , emit: tree
    tuple val(meta), path("*.iqtree")       , emit: iqtree_log
    tuple val(meta), path("*.log")          , emit: log
    tuple val(meta), path("*.mldist"), optional: true, emit: mldist
    tuple val(meta), path("*.model.gz"), optional: true, emit: model
    tuple val(meta), path("*.splits.nex"), optional: true, emit: splits
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${phylip.baseName}"
    def threads = task.ext.threads ?: task.cpus
    def sequence_type = task.ext.sequence_type ?: "DNA"
    // def model = task.ext.model ?: "TEST"
    def bootstrap = task.ext.bootstrap ?: "1000"
    def alrt = task.ext.alrt ?: "1000"
    def memory = task.ext.memory ? "-mem ${task.ext.memory}" : ""
    """
    iqtree \\
        -nt ${threads} \\
        -s ${phylip} \\
        -pre ${prefix} \\
        -st ${sequence_type} \\
        -bb ${bootstrap} \\
        -alrt ${alrt} \\
        ${memory} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree --version 2>&1 | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${phylip.baseName}"
    """
    touch ${prefix}.treefile
    touch ${prefix}.iqtree
    touch ${prefix}.log
    touch ${prefix}.mldist
    touch ${prefix}.model.gz
    touch ${prefix}.splits.nex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree --version 2>&1 | head -n1 | cut -d' ' -f3 || echo "2.3.4")
    END_VERSIONS
    """
}
