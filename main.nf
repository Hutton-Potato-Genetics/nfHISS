include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

help_message = """
Needin' some help? ------------------------------------------------------------

There are currently three workflows available in this pipeline:
- agrenseq: for carrying out association analysis from renseq data
- drenseq: for calculating read coverage of known resistance genes

agrenseq ----------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow agrenseq \
                                               --reference <reference> \
                                               --reads <read_scores> \
                                               <other options>

Options:
    --reference <reference>    Path to the reference fasta
    --reads <read_scores>      Path to the reads file - tab-separated file
                               with columns 'sample', 'forward', 'reverse',
                               and 'score'

drenseq -----------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow drenseq \
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

smrt-renseq --------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow smrtrenseq \
                                               --reads <reads_locations_tsv> \
                                               --genome_size
                                               <approximate_genome_size> \
                                               --max_input_coverage
                                               <maximum_input_coverage> \
                                               --flanking <flanking_bases>

Options:
    --reads <reads_locations_tsv>                   Path to the read locations
                                                    file - tab-separated file
                                                    with columns 'sample' and
                                                    'reads'
    --genome_size <approximate_genome_size>         Approximate expected
                                                    assembly size - parameter
                                                    for HiCanu
    --max_input_coverage <maximum_input_coverage    Maximum input coverage used
                                                    by HiCanu. If you have more
                                                    coverage than this, HiCanu
                                                    will randomly downsample
                                                    your input. Recommend
                                                    setting this to an excess
                                                    of your coverage.
    --flanking <flanking_bases>                     Number of bases to use as a
                                                    flanking region for NLR
                                                    Annotators fasta output
"""

workflow {
    switch (params.workflow) {
        case "agrenseq":
            agrenseq()
            break
        case "drenseq":
            drenseq()
            break
        case "smrtrenseq":
            smrtrenseq()
            break
        default:
            error("Unknown workflow: ${params.workflow}")
            break
    }
}
