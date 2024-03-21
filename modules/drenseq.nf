process Fastp {
    container 'https://depot.galaxyproject.org/singularity/fastp:0.23.3--h5f740d0_0'
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path('R1.fastq.gz'), path('R2.fastq.gz')
    script:
    """
    fastp -i $read1 -I $read2 -o R1.fastq.gz -O R2.fastq.gz
    """
}


process BowtieBuild {
    container 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h3321a2d_0'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '1h'
    input:
    path reference
    output:
    path 'bowtie2_index'
    script:
    """
    mkdir bowtie2_index
    bowtie2-build $reference bowtie2_index/reference
    """
}

process BowtieAlign {
    container 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h3321a2d_0'
    cpus 8
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    path bowtie2_index
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path('aligned.sam')
    script:
    """
    bowtie2 \
      -x ${bowtie2_index}/reference \
      -1 $read1 \
      -2 $read2 \
      --rg-id $sample \
      --rg SM:${sample} \
      -p ${task.cpus} \
      --score-min L,0,-0.24 \
      --phred33 \
      --fr \
      --maxins 1000 \
      --very-sensitive \
      --no-unal \
      --no-discordant \
      -k 10 \
      > aligned.sam
    """
}

process SamtoolsSort {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path("aligned.sorted.bam")
    script:
    """
    samtools sort -@ ${task.cpus} -o aligned.sorted.bam ${sam}
    """
}

process SambambaFilter {
    container 'https://depot.galaxyproject.org/singularity/sambamba:1.0--h98b6b92_0'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path('aligned.filtered.bam')
    script:
    """
    sambamba view --format=bam --filter='[NM] == 0' ${bam} > aligned.filtered.bam
    """
}

process BedtoolsCoverage {
    publishDir 'coverage', mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '1h'
    input:
    path bed
    tuple val(sample), path(bam)
    output:
    path '${sample}.coverage.txt'
    script:
    """
    bedtools coverage -a $bed -b $bam > ${sample}.coverage.txt
    """
}

process FreeBayes {
    container 'https://depot.galaxyproject.org/singularity/freebayes:1.3.7--h1870644_0'
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    path reference
    path bed
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path('variants.vcf')
    script:
    """
    freebayes \
      -f ${reference} \
      -t ${bed} \
      --min-alternate-count 2 \
      --min-alternate-fraction 0.05 \
      --ploidy 2 \
      --throw-away-indel-obs \
      --throw-away-mnps-obs \
      --throw-away-complex-obs \
      -m 0 \
      -v variants.vcf \
      --legacy-gls ${bam}
    """
}

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    bed = file(params.bed)

    reads = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) }

    trimmed_reads = Fastp(reads)
    bams = BowtieAlign(bowtie2_index.first(), trimmed_reads)
    sorted_bams = SamtoolsSort(bams)
    filtered_bams = SambambaFilter(sorted_bams)
    BedtoolsCoverage(bed, filtered_bams)
}
