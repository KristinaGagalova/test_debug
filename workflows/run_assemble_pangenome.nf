include { ASSEMBLE_SHORT_READS } from '../subworkflows/local/assemble_short_reads/main'
include { ASSEMBLE_NANO } from '../subworkflows/local/assemble_nano/main'
include { POLISH_GENOME } from '../subworkflows/local/polish_genome/main'

if (params.in_table) {in_table = params.in_table} else { exit 1, 'No reads specified!' }

def samples = []

new File(params.in_table).eachLine { line, lineNumber ->
    if (lineNumber > 0) { 
        def columns = line.split("\t").collect { it ?: null }
        def samplename = columns[0] 
        def reads1 = columns[1] != null ? file(columns[1]) : null
        def reads2 = columns[2] != null ? file(columns[2]) : null
        def longreads = columns[3] != null ? file(columns[3]) : null
        def contam = columns[4] != null ? columns[4] : "False"

        if (reads1 != null) {
            if (columns[1].split(",").size() > 1) {
                reads1 = columns[1].split(',').collect {file(it)}
            }
        }

        if (reads2 != null) {
            if (columns[2].split(",").size() > 1) {
                reads2 = columns[2].split(',').collect {file(it)}
            }
        }

        if (longreads != null) {
            if (columns[3].split(",").size() > 1) {
                longreads = columns[3].split(',').collect {file(it)}
            }
        }

        sample = tuple(samplename, reads1, reads2, longreads, contam)
        
        samples << sample
    }
}


ch_samples = Channel.from(samples).branch {
    only_paired_reads: it[1] != null && it[2] != null && it[3] == null
    only_single_reads: it[1] != null && it[2] == null && it[3] == null
    has_long_reads: it[3] != null
    other: true
}
.set { branched_samples }



workflow ASSEMBLE_PANGENOME_WORKFLOW{
    ASSEMBLE_SHORT_READS(branched_samples.only_paired_reads)


    ASSEMBLE_NANO(branched_samples.has_long_reads)

    // branched_samples.has_long_reads.join(ASSEMBLE_NANO.out.assembly)
    // .set {assembled_long_reads}

    // assembled_long_reads.branch {
    //     only_long_reads: it[1] == null && it[2] == null && it[3] != null
    //     paired_and_long_reads: it[1] != null && it[2] != null && it[3] != null
    // }
    // .set {branched_assembled_long_reads}

    // POLISH_GENOME(branched_assembled_long_reads.paired_and_long_reads)


    // only_paired_reads: it[1] != null && it[2] != null && it[3] == null
    // paired_and_long_reads: it[1] != null && it[2] != null && it[3] != null
    // only_single_reads: it[1] != null && it[2] == null && it[3] == null
    // single_and_long_reads: it[1] != null && it[2] == null && it[3] != null
    // only_long_reads: it[1] == null && it[2] == null && it[3] != null
}