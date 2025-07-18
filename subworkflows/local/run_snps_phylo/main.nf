include { BWA_INDEX }                    from '../../../modules/local/software/bwa/indexgenome.nf'
include { SAMTOOLS_INDEX_GENOME }        from '../../../modules/local/software/samtools/indexgenome.nf'
include { DICT_REF }                     from '../../../modules/local/software/gatk/index_genome.nf'
include { MARK_ILLUMINA_ADAPTERS }       from '../../../modules/local/software/gatk/mark_adapters.nf'
include { BWA_MAPREADS as ALIGN_TO_REF_UBAM } from '../../../modules/local/software/bwa/mapreads.nf'
include { INDEX_BAM }                    from '../../../modules/local/software/samtools/index_bam.nf'
include { MARK_DUPLICATES }              from '../../../modules/local/software/gatk/mark_duplicates.nf'
include { MERGE_BAM_WITH_UBAM }          from '../../../modules/local/software/gatk/merge_bam_with_ubam.nf'
include { COMBINE_AND_GENOTYPE_VCF }     from '../../../modules/local/software/gatk/combine_and_genotype_vcf.nf'
include { FILTER_SNPS_AND_INDELS }       from '../../../modules/local/software/gatk/filter_snps_and_indels.nf'
include { QUALITY_FILTER_VARIANTS }      from '../../../modules/local/software/gatk/quality_filter_variants.nf'
include { FINAL_FILTER_VARIANTS }        from '../../../modules/local/software/bcftools/final_filter_variants.nf'
include { VCF_TO_PHYLIP }                from '../../../modules/local/software/scripting/vcf2phylip.nf'
include { GENERATE_TREE }                from '../../../modules/local/software/iqtree/generate_tree.nf' 

workflow CREATE_INDEX {
    take: 
        reference // tuple(meta, fasta)

    main:
        BWA_INDEX(reference)
        SAMTOOLS_INDEX_GENOME(reference)
        DICT_REF(reference)

        // Combine all reference-related files
        ref_bundle = reference
            .join(SAMTOOLS_INDEX_GENOME.out.index, by: 0)
            .join(DICT_REF.out.index, by: 0)
            .map { meta, fasta, fai, dict -> 
                [meta, fasta, fai, dict] 
            }

    emit:
        bwa_index = BWA_INDEX.out.index      // tuple(meta, index)
        seq_dict = DICT_REF.out.index        // tuple(meta, index)
        ref_bundle = ref_bundle              // tuple(meta, fasta, fai, dict)
}

workflow UBAM_QC_AND_MAPPING {
    take:
        reads      // tuple(meta, reads)
        ubams      // tuple(meta, ubam)
        reference  // tuple(meta, reference)
        ref_bundle // tuple(meta, fasta, fai, dict)
        bwa_index  // tuple(meta, index)

    main:
        MARK_ILLUMINA_ADAPTERS(reads, ubams)
        ALIGN_TO_REF_UBAM(MARK_ILLUMINA_ADAPTERS.out.marked_fastq, bwa_index, reference)
        INDEX_BAM(ALIGN_TO_REF_UBAM.out.bam)
        MARK_DUPLICATES(INDEX_BAM.out.bam)

        // Properly pair the BAM files with their corresponding uBAM files by sample ID
        paired_bams = MARK_DUPLICATES.out.bam
            .join(ubams, by: 0)  // Join by meta (index 0)
            .map { meta, bam, ubam -> [meta, bam, ubam] }

        MERGE_BAM_WITH_UBAM(paired_bams, ref_bundle)
    
    emit:
        ubam      = MARK_ILLUMINA_ADAPTERS.out.marked_ubam      // tuple(meta, ubam)
        bam       = MERGE_BAM_WITH_UBAM.out.bam                 // tuple(meta, bam)
        bam_index = MERGE_BAM_WITH_UBAM.out.bai                 // tuple(meta, bai)
}

workflow VCF_GENOTYPING_AND_FILTERING {
    take:
        gvcfs     // tuple(meta, gvcf)
        reference // tuple(meta, fasta)
        fai_index // tuple(meta, index)
        // seq_dict  // tuple(meta, seq_dict)

    main:
        COMBINE_AND_GENOTYPE_VCF(gvcfs, reference, fai_index)
        FILTER_SNPS_AND_INDELS(COMBINE_AND_GENOTYPE_VCF.out.vcf, \
                                COMBINE_AND_GENOTYPE_VCF.out.vcf_index, 
                                reference, fai_index)
        QUALITY_FILTER_VARIANTS(FILTER_SNPS_AND_INDELS.out.snps_vcf,
                                FILTER_SNPS_AND_INDELS.out.snps_vcf_index,
                                FILTER_SNPS_AND_INDELS.out.indels_vcf,
                                FILTER_SNPS_AND_INDELS.out.indels_vcf_index,
                                reference,
                                fai_index)
        FINAL_FILTER_VARIANTS(QUALITY_FILTER_VARIANTS.out.filtered_snps_vcf,
			        QUALITY_FILTER_VARIANTS.out.filtered_snps_vcf_index,
				QUALITY_FILTER_VARIANTS.out.filtered_indels_vcf,
				QUALITY_FILTER_VARIANTS.out.filtered_indels_vcf_index)
    
    emit:
        vcf = FINAL_FILTER_VARIANTS.out.vcf // tuple(meta, vcf)
}

workflow BUILD_TREE {
    take:
        vcf // tuple(meta, vcf)

    main:
        VCF_TO_PHYLIP(vcf)
        GENERATE_TREE(VCF_TO_PHYLIP.out.fasta)
    
    emit:
        tree = GENERATE_TREE.out.tree // tuple(meta, tree)
}
