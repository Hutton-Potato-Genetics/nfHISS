include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    if (params.workflow == "agrenseq") {
        agrenseq()
    } else if (params.workflow == "drenseq") {
        drenseq()
    } else if (params.workflow == "smrtrenseq") {
        smrtrenseq()
    } else {
        error("Unknown workflow: ${params.workflow}")
    }
}
