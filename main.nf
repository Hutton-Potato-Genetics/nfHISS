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
                                --reads <read_scores> \
                                <other options>

Options:
  --reference <reference>     Path to the reference fasta
  --reads <read_scores>       Path to the reads file - tab-separated file
                              with columns 'sample', 'forward', 'reverse',
                              and 'score'

drenseq -----------------------------------------------------------------------

Usage:
  nextflow run SwiftSeal/nfHISS --workflow drenseq \
                                --reference <reference> \
                                --bed <bed_file> \
                                --reads <read_scores> \
                                <other options>

Options:
  --reference <reference>     Path to the reference fasta
  --reads <read_scores>       Path to the reads file - tab-separated file
  --bed <bed_file>            Path to the bed file - tab-separated file
                              with columns 'sample', 'forward', and
                              'reverse'
"""

workflow {
    switch (params.workflow) {
        case "agrenseq":
            agrenseq()
            break
        case "drenseq":
            drenseq()
            break
        default:
            error("Unknown workflow: ${params.workflow}")
            break
    }
}
