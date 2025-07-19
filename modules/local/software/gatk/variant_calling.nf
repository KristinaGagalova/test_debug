process VARIANT_CALLING {

    tag "${meta.id}"
    label 'very_large_task'
    label 'gatk'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'broadinstitute/gatk:4.2.6.1' :
        'broadinstitute/gatk:4.2.6.1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta_bai), path(bam_index)
    // ----
    tuple val(meta_ref), path(reference)
    tuple val(meta_fai), path(fai_index)
    tuple val(meta_dict), path(seq_dict)

    output:
    tuple val(meta), path("${meta.id}.g.vcf"), emit: gvcf
    tuple val(meta), path("${meta.id}.g.vcf.idx"), emit: gvcf_index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_mapping_quality = task.ext.min_mapping_quality ?: 20
    def min_base_quality = task.ext.min_base_quality ?: 20
    def native_pair_hmm_threads = task.ext.native_pair_hmm_threads ?: Math.max(1, task.cpus)
    
    """
    mkdir -p tmp

    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        --input ${bam} \\
        --output ${prefix}.g.vcf \\
        --reference ${reference} \\
        --emit-ref-confidence GVCF \\
        --minimum-mapping-quality ${min_mapping_quality} \\
        --min-base-quality-score ${min_base_quality} \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --annotation-group StandardHCAnnotation \\
        --tmp-dir tmp \\
        --native-pair-hmm-threads ${native_pair_hmm_threads} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -E '^(GATK|The Genome Analysis Toolkit)' | sed 's/.*v//' | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.g.vcf
    touch ${prefix}.g.vcf.idx
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}
