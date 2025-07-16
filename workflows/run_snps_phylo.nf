include { CREATE_INDEX }                     from '../subworkflows/local/run_snps_phylo/main.nf'
include { UBAM_QC_AND_MAPPING }              from '../subworkflows/local/run_snps_phylo/main.nf'
include { VCF_GENOTYPING_AND_FILTERING }     from '../subworkflows/local/run_snps_phylo/main.nf'
include { BUILD_TREE }                       from '../subworkflows/local/run_snps_phylo/main.nf'

include { GENERATE_UBAM }                    from '../modules/local/software/gatk/generate_ubam.nf'
include { VARIANT_CALLING }                  from '../modules/local/software/gatk/variant_calling.nf'

// Parameter validation
if (params.ref) { 
    reference = file(params.ref, checkIfExists: true) 
} else { 
    exit 1, 'No reference genome specified!' 
}

if (params.fastq) {
    fastq = Channel.fromFilePairs(params.fastq, checkIfExists: true)
        .map { sample_id, reads -> 
            def meta = [id: sample_id]
            [meta, reads]
        }
} else { 
    exit 1, 'No reads specified!' 
}

workflow SNPS_PHYLO_WORKFLOW {
    // Create reference metadata
    ref_meta = [id: reference.baseName]
    ref_tuple = tuple(ref_meta, reference)
    
    // Create index files
    CREATE_INDEX(ref_tuple)
    
    // Generate uBAM files
    GENERATE_UBAM(fastq)
    
    // QC and mapping
    UBAM_QC_AND_MAPPING(
        GENERATE_UBAM.out, 
        ref_tuple, 
        CREATE_INDEX.out.bwa_index, 
        CREATE_INDEX.out.seq_dict
    )
    
    // Variant calling
    VARIANT_CALLING(
        UBAM_QC_AND_MAPPING.out, 
        ref_tuple, 
        CREATE_INDEX.out.fai_index, 
        CREATE_INDEX.out.seq_dict
    )
    
    // Collect all variant calls and process
    VCF_GENOTYPING_AND_FILTERING(
        VARIANT_CALLING.out.collect(), 
        ref_tuple, 
        CREATE_INDEX.out.fai_index, 
        CREATE_INDEX.out.seq_dict
    )
    
    // Build phylogenetic tree
    BUILD_TREE(VCF_GENOTYPING_AND_FILTERING.out)
}