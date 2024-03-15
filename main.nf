include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'

help_message = """
Needin' some help? ------------------------------------------------------------

There are currently two workflows available in this pipeline:
- agrenseq: for carrying out association analysis from renseq data
- drenseq: for calculating read coverage of known resistance genes

agrenseq ----------------------------------------------------------------------

Usage:
    nextflow run SwiftSeal/nfHISS --workflow agrenseq \
                                  --reference <reference> \
                                  --read_scores <read_scores> \
                                  <other options>

Options:
    --reference <reference>     Path to the reference fasta
    --read_scores <read_scores> Path to the read scores file
"""

workflow {
    switch (params.workflow) {
        case "agrenseq":
            agrenseq()
        case "drenseq":
            drenseq()
        default:
            error("Unknown workflow: ${params.workflow}")
    }
}
