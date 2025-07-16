include { KRAKEN2_SETUP }           from '../../../modules/local/software/kraken/kraken2setup'
include { KRAKEN2 }                 from '../../../modules/local/software/kraken/kraken2'
include { EXTRACT_KRAKEN_READS }    from '../../../modules/local/software/kraken/krakentools'
include { BBSPLIT }                 from '../../../modules/local/software/bbmap/bbsplit'

workflow SHORT_READ_CONTAM_FILTERING {

    take:
    reads

    main:
    reference = null
    bin = null
    kraken_db = null

    if (params.ref) { reference = Channel.fromPath(params.ref).collect() } else { println 'No reference genome/s specified! (ignore if not performing contamination filtering)' }
    if (params.bin) { bin = params.bin } else { println 'No bin selection specified! (ignore if not performing contamination filtering)' }
    if (params.kraken_db) { kraken_db = params.kraken_db } else { println 'No kraken database specified! (ignore if not performing contamination filtering)' }

    if (reference && bin && kraken_db) {
        KRAKEN2(reads, kraken_db)
        .set{ kraken_results }

        reads.join(kraken_results.kraken_output)
        .set { reads_and_kraken }

        EXTRACT_KRAKEN_READS(reads_and_kraken)
        .set { extracted_reads }

        BBSPLIT(extracted_reads.cleaned_reads, reference, bin)
        .set { binned_reads }

        final_bin = binned_reads.reads_in_bin
        }
    
    else {
        final_bin = null
        }

    emit:
    contam_filtered_reads = final_bin
    
}