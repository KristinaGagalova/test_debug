include { BWA_INDEX }                    from '../../../modules/local/software/bwa/indexgenome.nf'
include { SAMTOOLS_INDEX_GENOME }        from '../../../modules/local/software/samtools/indexgenome.nf'
include { DICT_REF }                     from '../../../modules/local/software/gatk/index_genome.nf'
include { MARK_ILLUMINA_ADAPTERS }       from '../../../modules/local/software/gatk/mark_adapters.nf'
include { BWA_MAPREADS as ALIGN_TO_REF_UBAM } from '../../../modules/local/software/bwa/mapreads.nf'
include { SORT_AND_INDEX_BAM }           from '../../../modules/local/software/samtools/sort_and_index_bam.nf'
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

    emit:
        bwa_index = BWA_INDEX.out.index            // tuple(meta), path(index) path to files
        seq_dict = DICT_REF.out.index              // tuple(meta), path(index) path to index
        fa_index = SAMTOOLS_INDEX_GENOME.out.index // tuple(meta), path(index) path to index
}

workflow UBAM_QC_AND_MAPPING {
    take:
        reads      // tuple(meta, reads)
        ubams      // tuple(meta, ubam)
        reference  // path(reference)
        fa_index   // path(fai)
        bwa_index  // path(index)
        seq_dict   // path(seq_dict)

    main:
        MARK_ILLUMINA_ADAPTERS(reads, ubams)

        // Extract values from tuples - these are already channels
        bwa_index_val = bwa_index.map { meta, index -> index }.first()
        reference_val = reference.map { meta, fasta -> fasta }.first()
        fa_index_val = fa_index.map { meta, index -> index }.first()
        seq_dict_val = seq_dict.map { meta, index -> index }.first()

        ALIGN_TO_REF_UBAM(MARK_ILLUMINA_ADAPTERS.out.marked_fastq,
				bwa_index_val,
				reference_val)
        SORT_AND_INDEX_BAM(ALIGN_TO_REF_UBAM.out.bam)
        MARK_DUPLICATES(SORT_AND_INDEX_BAM.out.bam)

        // Pair the BAM files with their corresponding uBAM files by sample ID
        paired_bams = MARK_DUPLICATES.out.bam
            .join(MARK_ILLUMINA_ADAPTERS.out.marked_ubam, by: 0)  // Join by meta (index 0)
	    .join(MARK_DUPLICATES.out.index, by: 0)  // Join the BAM index file
            .map { meta, bam, ubam, bai -> [meta, bam, ubam, bai] }

        MERGE_BAM_WITH_UBAM(paired_bams, reference_val, fa_index_val, seq_dict_val)
    
    emit:
        ubam      = MARK_ILLUMINA_ADAPTERS.out.marked_ubam      // tuple(meta, ubam)
        bam       = MERGE_BAM_WITH_UBAM.out.bam                 // tuple(meta, bam, bai)
        //bam_index = MERGE_BAM_WITH_UBAM.out.bai                 // tuple(meta, bai)
}

workflow VCF_GENOTYPING_AND_FILTERING {
    take:
        meta      // meta value for sample name
        gvcfs     // tuple(meta, gvcf) # collected gvcfs 
        reference // tuple(meta, fasta)
        fa_index  // tuple(meta, index)
        seq_dict  // tuple(meta, seq_dict)

    main:

        // Extract values from tuples - these are already channels
        reference_val = reference.map { meta, fasta -> fasta }.first()
        fa_index_val = fa_index.map { meta, index -> index }.first()
        seq_dict_val = seq_dict.map { meta, index -> index }.first()
        
        COMBINE_AND_GENOTYPE_VCF(meta, gvcfs, reference_val, fa_index_val, seq_dict_val)
        // reintroduce meta id
        transformed_channel = COMBINE_AND_GENOTYPE_VCF.out.vcf.map { meta_string, vcf, idx ->
 		   	def meta_map = [id: meta_string]
    		   	[meta_map, vcf, idx]
			}
        FILTER_SNPS_AND_INDELS(transformed_channel,
                                reference_val,
				fa_index_val,
				seq_dict_val)
        QUALITY_FILTER_VARIANTS(FILTER_SNPS_AND_INDELS.out.vcfs,
                                reference_val,
                                fa_index_val,
				seq_dict_val)
        FINAL_FILTER_VARIANTS(QUALITY_FILTER_VARIANTS.out.filtered_vcfs)
    
    emit:
        vcf = FINAL_FILTER_VARIANTS.out.vcf // tuple(meta, vcf)
}

workflow BUILD_TREE {
    take:
        vcf // tuple(meta, vcf)

    main:
        VCF_TO_PHYLIP(vcf)
        GENERATE_TREE(VCF_TO_PHYLIP.out.phylip)
    
    emit:
        tree = GENERATE_TREE.out.tree // tuple(meta, tree)
}
