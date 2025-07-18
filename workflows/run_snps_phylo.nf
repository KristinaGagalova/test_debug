include { CREATE_INDEX }                     from '../subworkflows/local/run_snps_phylo/main.nf'
include { UBAM_QC_AND_MAPPING }              from '../subworkflows/local/run_snps_phylo/main.nf'
include { VCF_GENOTYPING_AND_FILTERING }     from '../subworkflows/local/run_snps_phylo/main.nf'
include { BUILD_TREE }                       from '../subworkflows/local/run_snps_phylo/main.nf'

include { GENERATE_UBAM }                    from '../modules/local/software/gatk/generate_ubam.nf'
include { VARIANT_CALLING }                  from '../modules/local/software/gatk/variant_calling.nf'

workflow SNPS_PHYLO_WORKFLOW {

    // Parameter validation and input channels
    if (!params.ref) {
        exit 1, 'No reference genome specified!'
    }
    if (!params.fastq) {
        exit 1, 'No reads specified!'
    }

    reference = file(params.ref, checkIfExists: true)
    fastq = Channel.fromFilePairs(params.fastq, checkIfExists: true)
        .map { sample_id, reads -> 
            def meta = [id: sample_id]
            [meta, reads]
        }

    // Create reference metadata
    ref_tuple = Channel.of([
                [id: reference.baseName], 
                reference
                ])
    
    // Create index files
    CREATE_INDEX(ref_tuple)
    
    // Generate uBAM files
    GENERATE_UBAM(fastq)

    // QC and mapping
    UBAM_QC_AND_MAPPING(
        fastq,
        GENERATE_UBAM.out.ubam, 
        ref_tuple, 
        CREATE_INDEX.out.ref_bundle,
        CREATE_INDEX.out.bwa_index
    )
    
        
    // Variant calling
    //VARIANT_CALLING(
    //    UBAM_QC_AND_MAPPING.out.bam,
    //    UBAM_QC_AND_MAPPING.out.bam_index,
    //    ref_tuple
    //)
    
    // Collect all variant calls and process
    //VCF_GENOTYPING_AND_FILTERING(
    //    VARIANT_CALLING.out.gvcf.collect(), 
    //    ref_tuple, 
    //    CREATE_INDEX.out.fai_index
    //)
    
    // Build phylogenetic tree
    //BUILD_TREE(VCF_GENOTYPING_AND_FILTERING.out.vcf)
}
