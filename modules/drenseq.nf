process BowtieBuild {
    container 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h3321a2d_0'
    cpus 1
    memory { 1.GB * task.attempt }
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
    time '4h'
    input:
    path bowtie2_index
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path("${sample}.sam")
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
      > ${sample}.sam
    """
}

process SamtoolsSort {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0'
    cpus 1
    memory { 1.GB * task.attempt }
    time '4h'
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path("${sample}.sorted.bam")
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sam}
    """
}

process SambambaFilter {
    container 'https://depot.galaxyproject.org/singularity/sambamba:1.0--h98b6b92_0'
    cpus 1
    memory { 1.GB * task.attempt }
    time '4h'
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path("${sample}.filtered.bam")
    script:
    """
    sambamba view --format=bam --filter='[NM] == 0' ${bam} > ${sample}.filtered.bam
    """
}

process BedtoolsCoverage {
    publishDir 'coverage', mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3'
    cpus 1
    memory { 1.GB * task.attempt }
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

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    bed = file(params.bed)

    Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) }
        .set { reads }

    bams = BowtieAlign(bowtie2_index.first(), reads) \
        | SamtoolsSort \
        | SambambaFilter

    BedtoolsCoverage(bed, bams)
}
