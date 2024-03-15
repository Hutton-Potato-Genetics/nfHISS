process JellyfishCount {
    conda 'jellyfish=2.2.10'
    cpus 2
    memory '16 GB'
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path("${sample}.count")
    script:
    """
    zcat $read1 $read2 | jellyfish count -C -m 51 -s 1G -t ${task.cpus} -o ${sample}.count /dev/fd/0
    """
}

process JellyfishDump {
    conda 'jellyfish=2.2.10'
    cpus 1
    memory '2 GB'
    input:
    tuple val(sample), path(input)
    output:
    tuple val(sample), path("${sample}.dump")
    script:
    """
    jellyfish dump -L 10 -ct $input > ${sample}.dump
    """
}

process CreatePresenceMatrix {
    cpus 1
    memory '32 GB'
    input:
    path accession_table
    output:
    path 'presence_matrix.txt'
    script:
    """
    CreatePresenceMatrix.sh -i $accession_table -o presence_matrix.txt
    """
}

process NLRParser {
    conda 'meme=5.4.1=py310pl5321h9f004f7_2'
    cpus 4
    memory '4 GB'
    input:
    path assembly
    output:
    path 'nlrparser.txt'
    script:
    """
    NLR_Parser.sh -t ${task.cpus} -i $assembly -o nlrparser.txt
    """
}

process RunAssociation {
    cpus 1
    memory '16 GB'
    input:
    path presence_matrix
    path assembly
    path phenotype
    path nlrparser
    output:
    path 'agrenseq_result.txt'
    script:
    """
    RunAssociation.sh -i $presence_matrix -n $nlrparser -p $phenotype -a $assembly -o agrenseq_result.txt
    """
}

process BowtieBuild {
    container 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h3321a2d_0'
    cpus 1
    memory '2 GB'
    time '4h'
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
    time '4h'
    memory '8 GB'
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
    memory '4 GB'
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
    memory '4 GB'
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
    container 'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3'
    cpus 1
    memory '2 GB'
    time '1h'
    input:
    path bed
    tuple val(sample), path(bam)
    output:
    stdout
    script:
    """
    bedtools coverage -a $bed -b $bam | \
      awk 'BEGIN {OFS="\\t"} {print "${sample}", \$0}'
    """
}

workflow agrenseq {
    accession_table = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) } \
        | JellyfishCount \
        | JellyfishDump \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(name: "accession.tsv", newLine: true)

    phenotype_file = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> "${row.sample}\t${row.score}" } \
        | collectFile(name: "phenotype.tsv", newLine: true)

    matrix = CreatePresenceMatrix(accession_table)

    assembly = Channel
        .fromPath(params.assembly)

    nlrparser = NLRParser(assembly)

    association = RunAssociation(matrix, assembly, phenotype_file, nlrparser)
}

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    bed = Channel.fromPath(params.bed)

    Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) }
        .set { reads }

    bams = BowtieAlign(bowtie2_index.first(), reads) \
        | SamtoolsSort \
        | SambambaFilter

    BedtoolsCoverage(bed.first(), bams).collectFile(name: 'coverage.txt')
}

workflow {
    switch (params.workflow) {
        case "agrenseq": agrenseq()
        case "drenseq": drenseq()
        default: fail("Unknown workflow: ${params.workflow}")
    }
}
