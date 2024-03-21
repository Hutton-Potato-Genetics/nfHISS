# nfHISS

nfHISS is a re-implementation of the HISS pipeline using Nextflow.
This has been created as a result of recent changes to Snakemake which have reduced its compatibility with SLURM.

## Running nfHISS

To run nfHISS, you will first need to have Nextflow installed.

All nfHISS pipelines are executed through a single command:

```
nextflow run SwiftSeal/nfHISS --workflow <workflow> <additional arguments>
```

This will download the latest version of nfHISS.
