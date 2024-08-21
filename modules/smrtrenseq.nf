process TrimReads {
    container 'docker://quay.io/biocontainers/cutadapt:4.9--py312hf67a6ed_0'
    scratch true
    cpus 8
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(reads)
    val five_prime
    val three_prime
    output:
    path "${sample}_trimmed.fastq.gz"
    script
    """
    cutadapt -j 8 -g ^$five_prime -a $three_prime\$ -o ${sample}_trimmed.fastq.gz $reads
    """
}

process Canu {
    container 'docker://swiftseal/smrtrenseq:latest'
    cpus 8
    memory { 16.GB * task.attempt }
    maxRetries 3
    time '24h'
    input:
    tuple val(sample), path(reads)
    output:
    tuple val(sample), path("assembly/assembly.contigs.fasta")
    script:
    """
    canu \
        -d assembly \
        -p assembly \
        genomeSize=70m
        useGrid=false \
        -pacbio-hifi $reads \
        maxInputCoverage=20000
    """
}

process SeqkitStats {
    container 'docker://swiftseal/smrtrenseq:latest'
    cpus 1
    memory '1 GB'
    maxRetries 3
    time '1h'
    input:
    path reads
    output:
    path "${reads}.stats"
    script:
    """
    seqkit stats $reads > ${reads}.stats
    """
}

workflow smrtrenseq {
    trimmed_reads = Cutadapt(reads)
    assembly = Canu(trimmed_reads)
}
