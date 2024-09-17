include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

help_message = """
Needin' some help? ------------------------------------------------------------

There are currently three workflows available in this pipeline:
- smrtrenseq: for assembling HiFi RenSeq reads into contigs
- agrenseq: for carrying out association analysis from renseq data
- drenseq: for calculating read coverage of known resistance genes

agrenseq ----------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow agrenseq \
                                               --association_reference <association_reference> \
                                               --reads <read_scores> \
                                               --adaptor_1 <barcode_fasta_1> \
                                               --adaptor_2 <barcode_fasta_2> \
                                               --blast_reference <blast_reference> \
                                               --threshold <association_threshold> \
                                               --title <plot_title>



Options:
    --association_reference <association_reference>     Path to the association
                                                        reference fasta
    --reads <read_scores>                               Path to the reads file 
                                                        - tab-separated file
                                                        with columns 'sample'
                                                        'forward', 'reverse',
                                                        and 'score'
    --adaptor_1 <barcode_fasta_1>                       Path to first barcode
                                                        fasta
    --adaptor_2 <barcode_fasta_2>                       Path to second barcode
                                                        fasta
    --blast_reference <blast_reference>                 Path to the blast
                                                        reference fasta
    --threshold <association_threshold>                 Value to use as
                                                        significance threshold
                                                        in the association plot
    --title <plot_title>                                Title for association
                                                        plot
    --annotator_bed <bed_file_of_nlrs>                  Bed file of NLR
                                                        locations used to
                                                        identify candidate NLRs

drenseq -----------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow drenseq \
                                               --reference <reference> \
                                               --bed <bed_file> \
                                               --reads <read_scores> \
                                               <other options>

Options:
    --reference <reference>                     Path to the reference fasta
    --reads <read_scores>                       Path to the reads file -
                                                tab-separated file
    --bed <bed_file>                            Path to the bed file -
                                                tab-separated file with columns
                                                'sample', 'forward', and
                                                'reverse'
    --adaptor_1 <barcode_fasta_1>               Path to first barcode fasta
    --adaptor_2 <barcode_fasta_2>               Path to second barcode fasta
    --score <bowtie2_score_min>                 Parameter for BowTie2 to
                                                control allowed mismatch rate
    --max_align <maximum_allowed_alignments>    Parameter for BowTie2 to control
                                                maximum alignments allowed
    --baits <renseq_baits>                      Path to fasta file of RenSeq
                                                baits
    --identity <percent_identity>               Parameter for blastn - minimum
                                                percentage identity
    --coverage <coverage>                       Parameter for blastn - minimum
                                                coverage of hit
    --flank <flanking_region>                   Number of bases to take either
                                                side of a blastn hit
    --ulimit <ulimit>                           

smrt-renseq --------------------------------------------------------------------

Usage:
    nextflow run Hutton-Potato-Genetics/nfHISS --workflow smrtrenseq \
                                               --reads <reads_locations_tsv> \
                                               --genome_size
                                               <approximate_genome_size> \
                                               --max_input_coverage
                                               <maximum_input_coverage> \
                                               --flanking <flanking_bases> \
                                               --five_prime <5'_to_trim> \
                                               --three_prime <3'_to_trim>

Options:
    --reads <reads_locations_tsv>                   Path to the read locations
                                                    file - tab-separated file
                                                    with columns 'sample' and
                                                    'reads'
    --genome_size <approximate_genome_size>         Approximate expected
                                                    assembly size - parameter
                                                    for HiCanu
    --max_input_coverage <maximum_input_coverage>   Maximum input coverage used
                                                    by HiCanu. If you have more
                                                    coverage than this, HiCanu
                                                    will randomly downsample
                                                    your input. Recommend
                                                    setting this to an excess
                                                    of your coverage.
    --flanking <flanking_bases>                     Number of bases to use as a
                                                    flanking region for NLR
                                                    Annotators fasta output
    --five_prime <5'_to_trim>                       Sequence to be trimmed from
                                                    5' end of reads
    --three_prime <3'_to_trim>                      Sequence to be trimmed from
                                                    3' end of reads
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
