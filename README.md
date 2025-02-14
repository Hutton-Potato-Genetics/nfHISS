<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/HISSlogo_light.png">
  <img alt="Logo" src="assets/HISSlogo_dark.png">
</picture>

[![DOI:10.1186/s12859-023-05335-88](http://img.shields.io/badge/DOI-10.1186/s12859.023.05335.8-B31B1b.svg)](https://doi.org/10.1186/s12859-023-05335-8)
[![DOI](https://zenodo.org/badge/801906921.svg)](https://zenodo.org/doi/10.5281/zenodo.13789522)

# nfHISS

nfHISS is a re-implementation of the [HISS pipeline](https://github.com/SwiftSeal/HISS) using Nextflow.
This has been created as a result of recent changes to Snakemake which have reduced its compatibility with SLURM. Additionally a change has been made to favour Apptainer over Conda due to reported performance issues and some difficult to reproduce errors during enviornment resolution.

## Running nfHISS

To run nfHISS, you will first need to have [Nextflow installed](https://www.nextflow.io/docs/latest/install.html). Nextflow is also available on [bioconda](https://anaconda.org/bioconda/nextflow) for systems where users do not have sudo rights.

All nfHISS pipelines are executed through a single command:

```
nextflow run Hutton-Potato-Genetics/nfHISS -r main --workflow <workflow> <additional arguments>
```

This will download the latest version of nfHISS.

An example run through of the pipeline is provided, please use this if you intend to contribute or before reporting any errors. This also provides parameters used for tetraploid potato.

## Detailed options

All paths MUST be absolute paths.

### smrtrenseq

```
--reads <reads_locations_tsv>                   Path to read locations file - tab separated file
                                                    with columns 'sample' and 'reads'
--genome_size <appoximate_genome_size>          Approximate assembly size passed to HiCanu
--max_input_coverage <maximum input coverage>   Maximum input coverage passed to HiCanu
--flanking <flanking_bases>                     Flanking region for producing a fasta file with
                                                    NLR Annotator
--five_prime <5'_to_trim>                       5' sequence to be trimmed with Cutadapt
--three_prime <3'_to_trim>                      3' sequence to be trimmed with Cutadapt
```

### agrenseq

```
--association_reference <association_reference>     Path to fasta file of assembled contigs for
                                                    association
--reads <read_scores>                               Path to the reads file - tab-separated file
                                                    with columns 'sample', 'R1, 'R2' and
                                                    'score'. Where score represents the
                                                    phenotype, negative is susceptible and
                                                    positive is resistant
--adaptor_1 <barcode_fasta_1>                       Path to first barcode fasta passed to
                                                        Cutadapt
--adaptor_2 <barcode_fasta_2>                       Path to second barcode fasta passed to
                                                        Cutadapt
--blast_reference <blast_reference>                 Path to fasta file used as a reference for
                                                        blastn
--threshold <association_threshold>                 Value to use as significance threshold for
                                                        the association plot
--title <plot_title>                                Title for the assocition plot
--annotator_bed <bed_file_of_nlrs>                  Path to bed file used to location candidate
                                                        NLRs
```

### drenseq

```
--reference <reference_fasta>               Path to reference fasta file of candidate sequences
--reads <reads_file>                        Path to reads file - tab-separated with columns
                                                'sample','FRead', 'RRead'
--bed <bed_file>                            Path to the bed file of candidate NLRs
--adaptor_1  <barcode_fasta_1>              Path to first barcode fasta passed to Cutadapt
--adaptor_2  <barcode_fasta_2>              Path to second barcode fasta passed to Cutadapt
--score <bowtie2_score_min>                 Parameter passed to bowtie2 to control mismatch rate
--max_align <maximum_allowed_alignments>    Parameter passed to bowtie2 to control maximum
                                                number of allowed alignments per read
--baits <renseq_baits>                      Path to fasta file of the sequences of renseq baits
--identity <percentage_identity>            Parameter passed to blastn to filter results by
                                                percentage identity
--coverage <coverage>                       Parameter passed to blastn to filter results by
                                                coverage
--flank <flanking>                          Number of bases to take flanking a blastn hit of
                                                baits
--ulimit <ulimit>                           Controls a core unix parameter of maximum file open
                                                limit to handle large runs
```
