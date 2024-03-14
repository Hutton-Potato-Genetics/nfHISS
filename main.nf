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

workflow {
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
