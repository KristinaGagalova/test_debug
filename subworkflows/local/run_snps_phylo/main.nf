include { BWA_INDEX } from '../../../modules/local/software/bwa/indexgenome.nf'
include { SAMTOOLS_INDEX_GENOME } from '../../../modules/local/software/samtools/indexgenome.nf'
include { DICT_REF } from '../../../modules/local/software/gatk/index_genome.nf'
include { MARK_ILLUMINA_ADAPTERS } from '../../../modules/local/software/gatk/mark_adapters.nf'
include { ALIGN_TO_REF_UBAM } from '../../../modules/local/software/bwa/mapreads.nf'
include { SORT_AND_INDEX_BAM } from '../../../modules/local/software/samtools/sort_and_index_bam.nf'
include { MARK_DUPLICATES } from '../../../modules/local/software/gatk/mark_duplicates.nf'
include { MERGE_BAM_WITH_UBAM } from '../../../modules/local/software/gatk/merge_bam_with_ubam.nf'
include { COMBINE_AND_GENOTYPE_VCF } from '../../../modules/local/software/gatk/combine_and_genotype_vcf.nf'
include { FILTER_SNPS_AND_INDELS } from '../../../modules/local/software/gatk/filter_snps_and_indels.nf'
include { QUALITY_FILTER_VARIANTS } from '../../../modules/local/software/gatk/quality_filter_variants.nf'
include { FINAL_FILTER_VARIANTS } from '../../../modules/local/software/bcftools/final_filter_variants.nf'
include { VCF_TO_PHYLIP } from '../../../modules/local/software/scripting/vcf2phylip.nf'
include { GENERATE_TREE } from '../../../modules/local/software/iqtree/generate_tree.nf' 



workflow CREATE_INDEX {
    take: 
        reference //tuple(meta, fasta)

    main:
        BWA_INDEX(reference)
        SAMTOOLS_INDEX_GENOME(reference)
        DICT_REF(reference)

    emit:
        bwa_index = BWA_INDEX.out.index
        fai_index = SAMTOOLS_INDEX_GENOME.out.index
        seq_dict = DICT_REF.out.index
}


workflow UBAM_QC_AND_MAPPING {
    take:
        reads //tuple(sampleID, reads)
        ubams
        reference //tuple(meta, fasta)
        bwa_index //tuple(meta, index)
        seq_dict //tuple(meta, index)

    main:
        MARK_ILLUMINA_ADAPTERS(reads, ubams)
        ALIGN_TO_REF_UBAM(MARK_ILLUMINA_ADAPTERS.out, reference, bwa_index)
        SORT_AND_INDEX_BAM(ALIGN_TO_REF_UBAM.out)
        MARK_DUPLICATES(SORT_AND_INDEX_BAM.out)
        MERGE_BAM_WITH_UBAM(MARK_DUPLICATES.out, reference, seq_dict)
    
    emit:
        ubam = MERGE_BAM_WITH_UBAM.out[0]
        bam = MERGE_BAM_WITH_UBAM.out[1]
        bam_index = MERGE_BAM_WITH_UBAM.out[2]   
}


workflow VCF_GENOTYPING_AND_FILTERING {
    take:
        gvcfs
        reference //tuple(meta, fasta)
        fai_index //tuple(meta, index)
        seq_dict //tuple(meta, index)

    main:
        COMBINE_AND_GENOTYPE_VCF(gvcfs, reference, fai_index, seq_dict)
        FILTER_SNPS_AND_INDELS(COMBINE_AND_GENOTYPE_VCF.out, reference, fai_index, seq_dict)
        QUALITY_FILTER_VARIANTS(FILTER_SNPS_AND_INDELS.out, reference, fai_index, seq_dict)
        FINAL_FILTER_VARIANTS(QUALITY_FILTER_VARIANTS.out)
    
    emit:
        FINAL_FILTER_VARIANTS.out
}


workflow BUILD_TREE {
    take:
        vcf

    main:
        VCF_TO_PHYLIP(vcf)
        GENERATE_TREE(VCF_TO_PHYLIP.out)
    
    emit:
        tree = GENERATE_TREE.out
}