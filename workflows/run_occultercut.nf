////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

include { OCCULTERCUT_RUN       } from '../subworkflows/local/run_occultercut/main'

workflow OCCULTERCUT_WORKFLOW {

        samples_ch = Channel.fromPath(params.inputdir, type: 'file', checkIfExists: true)
			.map { path ->
			def name = "${path.baseName}"
			tuple(name, path)
			}

	OCCULTERCUT_RUN(samples_ch)
        println("Occultercut workflow - executed")
}
