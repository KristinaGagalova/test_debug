process FINAL_FILTER_VARIANTS {

    tag "${meta.id}"
    label 'small_task'
    label 'bcftools'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'staphb/bcftools:1.17':
        null }"

    input:
    tuple val(meta), path(snps)
    tuple val(meta), path(snps_idx)
    tuple val(meta), path(indels)
    tuple val(meta), path(indels_idx)

    output:
    tuple val(meta), path("*.filtered.biallelic.*.snps.vcf"), emit: vcf
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''  // bcftools view args
    def args2 = task.ext.args2 ?: ''  // bcftools +prune args
    def ref_name = params.refname ?: meta.id
    def missing_data = params.missingdata ?: task.ext.missing_data ?: "0.5"
    def max_alleles = task.ext.max_alleles ?: "2"
    def ld_threshold = task.ext.ld_threshold ?: "0.9"
    def window_size = task.ext.window_size ?: "5000bp"
    def output_format = task.ext.output_format ?: "v"
    
    // Create filename components with actual values
    def missing_suffix = "miss${missing_data.toString().replace('.', '')}"
    def window_suffix = window_size.replace("bp", "").replace("kb", "k")
    def ld_suffix = "ld${ld_threshold.toString().replace('.', '')}"
    def alleles_suffix = "max${max_alleles}alleles"
    
    def output_filename = "${ref_name}.combined_panel.filtered.biallelic.${alleles_suffix}.${missing_suffix}.1random${window_suffix}.${ld_suffix}.snps.vcf"
    """
    bcftools view \\
        --max-alleles ${max_alleles} \\
        --exclude-types indels \\
        -e "F_MISSING>=${missing_data}||FILTER!='PASS'" \\
        ${args1} \\
        ${snps} \\
        -o ${ref_name}.combined_panel.filtered.biallelic.snps.bcf
	
    bcftools +prune \\
        -m ${ld_threshold} \\
        -w ${window_size} \\
        ${ref_name}.combined_panel.filtered.biallelic.snps.bcf \\
        -O${output_format} \\
        -o ${output_filename} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def ref_name = params.refname ?: meta.id
    def missing_data = params.missingdata ?: task.ext.missing_data ?: "0.5"
    def max_alleles = task.ext.max_alleles ?: "2"
    def ld_threshold = task.ext.ld_threshold ?: "0.9"
    def window_size = task.ext.window_size ?: "5000bp"
    
    def missing_suffix = "miss${missing_data.toString().replace('.', '')}"
    def window_suffix = window_size.replace("bp", "").replace("kb", "k")
    def ld_suffix = "ld${ld_threshold.toString().replace('.', '')}"
    def alleles_suffix = "max${max_alleles}alleles"
    
    def output_filename = "${ref_name}.combined_panel.filtered.biallelic.${alleles_suffix}.${missing_suffix}.1random${window_suffix}.${ld_suffix}.snps.vcf"
    """
    touch ${output_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | cut -d' ' -f2 || echo "1.17")
    END_VERSIONS
    """
}