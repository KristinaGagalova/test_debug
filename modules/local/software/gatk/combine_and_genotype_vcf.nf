process COMBINE_AND_GENOTYPE_VCF {

    tag "${meta.id}"
    label 'mediumHighMem_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1':
        null }"

    input:
    path gvcfs
    tuple val(meta), path(reference)
    tuple val(meta1), path(fai_index)

    output:
    tuple val(meta), path("*.combined_panel.vcf")     , emit: vcf
    tuple val(meta), path("*.combined_panel.vcf.idx") , emit: vcf_index
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // CombineGVCFs args
    def args2 = task.ext.args2 ?: ''  // GenotypeGVCFs args
    def tmp_dir = task.ext.tmp_dir ?: "tmp"
    def ref_name = params.refname ?: meta.id
    """
    mkdir -p ${tmp_dir}
    
    ls *.g.vcf > list.gvcfs.list
    
    gatk CombineGVCFs \\
        -R ${reference} \\
        --variant list.gvcfs.list \\
        -O ${ref_name}.combined_panel.g.vcf \\
        --TMP_DIR ${tmp_dir} \\
        ${args1}
    
    gatk GenotypeGVCFs \\
        -R ${reference} \\
        --variant ${ref_name}.combined_panel.g.vcf \\
        -O ${ref_name}.combined_panel.vcf \\
        --TMP_DIR ${tmp_dir} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}')
    END_VERSIONS
    """

    stub:
    def ref_name = params.refname ?: meta.id
    """
    touch ${ref_name}.combined_panel.vcf
    touch ${ref_name}.combined_panel.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^The Genome Analysis Toolkit' | awk '{print \$6}' || echo "4.2.6.1")
    END_VERSIONS
    """
}