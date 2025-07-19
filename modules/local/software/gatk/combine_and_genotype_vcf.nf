process COMBINE_AND_GENOTYPE_VCF {

    tag "combined_panel"
    label 'mediumHighMem_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    val meta
    path gvcfs
    path(reference)
    path(fai_index)
    path(seq_dict)

    output:
    tuple val(meta), path("${meta}.vcf"), path("${meta}.vcf.idx"), emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // CombineGVCFs args
    def args2 = task.ext.args2 ?: ''  // GenotypeGVCFs args
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def prefix = "${meta}"
    """
    ls *.g.vcf > list.gvcfs.list
    
    gatk CombineGVCFs \\
        -R ${reference} \\
        --variant list.gvcfs.list \\
        -O ${prefix}.g.vcf \\
        ${args1}
    
    gatk GenotypeGVCFs \\
        -R ${reference} \\
        --variant ${prefix}.g.vcf \\
        -O ${prefix}.vcf \\
        --create-output-variant-index true \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """
}
