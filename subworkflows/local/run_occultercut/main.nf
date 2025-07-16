/*
 * Run Occultercut
 */

params.occultercut_options  = [:]
fasta = params.inputdir

include { OCCULTERCUT           } from '../../../modules/local/software/occultercut/main'
include { getOcculterCutRegionFrequencies } from '../../../modules/local/software/occultercut/auxiliary'
include { getOcculterCutGroupedRegionFrequencies } from '../../../modules/local/software/occultercut/auxiliary'
include { tidyOcculterCutGFFs } from '../../../modules/local/software/occultercut/auxiliary'

workflow OCCULTERCUT_RUN {

    take:
    fasta         // channel: tuple [ val(name), file(fasta) ] ]

    main:
    occulterCutRegions                       = Channel.empty()
    occulterCutGroupedRegions                = Channel.empty()
    occulterCutRegionFrequencies             = Channel.empty()
    getOcculterCutGroupedRegionFrequencies   = Channel.empty()
        
    OCCULTERCUT ( fasta ).occulterCutRegions.set { occulterCutRegions }
    occulterCutGroupedRegions     = OCCULTERCUT.out.occulterCutGroupedRegions
    
    occulterCutRegionsGen = occulterCutRegions.combine(fasta, by: 0)
    occulterCutGroupedRegionsGen = occulterCutGroupedRegions.combine(fasta, by: 0)

    getOcculterCutRegionFrequencies ( occulterCutRegionsGen )
							.occulterCutRegionFrequencies
							.set { occulterCutRegionFrequencies }

    getOcculterCutGroupedRegionFrequencies ( occulterCutGroupedRegionsGen )
                                                        .getOcculterCutGroupedRegionFrequencies
                                                        .set { occulterCutGroupedRegionFrequencies }
    
    occulterCutRegionFrequenciesMix = occulterCutRegionFrequencies.mix(occulterCutGroupedRegionFrequencies)
    tidyOcculterCutGFFs ( occulterCutRegionFrequenciesMix )
							.occulterCutGff
							.set { occulterCutTidy }
    
    emit:
    occulterCutRegions                // channel: [val(name), file("${name}_occultercut_grouped_regions.gff3")]
    occulterCutGroupedRegions         //channel: [val(name), file("${name}_occultercut_grouped_regions.gff3")]
    occulterCutTidy
}

