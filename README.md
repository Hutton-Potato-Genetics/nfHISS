<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/HISSlogo_light.png">
  <img alt="Logo" src="assets/HISSlogo_dark.png">
</picture>

[![DOI:10.1186/s12859-023-05335-88](http://img.shields.io/badge/DOI-10.1186/s12859.023.05335.8-B31B1b.svg)](https://doi.org/10.1186/s12859-023-05335-8)

# nfHISS

nfHISS is a re-implementation of the [HISS pipeline](https://github.com/SwiftSeal/HISS) using Nextflow.
This has been created as a result of recent changes to Snakemake which have reduced its compatibility with SLURM. Additionally a change has been made to favour Apptainer over Conda due to reported performance issues and some difficult to reproduce errors during enviornment resolution.

## Running nfHISS

To run nfHISS, you will first need to have [Nextflow installed](https://www.nextflow.io/docs/latest/install.html). Nextflow is also available on [bioconda](https://anaconda.org/bioconda/nextflow) for systems where users do not have sudo rights.

All nfHISS pipelines are executed through a single command:

```
nextflow run SwiftSeal/nfHISS --workflow <workflow> <additional arguments>
```

This will download the latest version of nfHISS.

An example run through of the pipeline is provided, please use this if you intend to contribute or before reporting any errors. This also provides parameters used for tetraploid potato.

## Detailed options

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
