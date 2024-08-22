process CountKmers {
    container 'docker://quay.io/biocontainers/kmc3.2.4--h6dccd9a_1'
    scratch true
    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStats == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path("${sample}.dump")
    script:
    """
    cat $read1 $read2 > reads.fq.gz
    kmc -k51 -m8 -t${task.cpus} reads.fq.gz kmc_output .
    kmc_tools transform kmc_output -ci10 dump ${sample}.dump
    """
}

process CreatePresenceMatrix {
    scratch true
    cpus 1
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStats == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
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
    container 'docker://quay.io/biocontainers/meme:5.5.6--pl5321h4242488_0'
    scratch true
    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStats == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '8h'
    input:
    path reference
    output:
    path 'nlrparser.txt'
    script:
    """
    nlr_parser.sh -t ${task.cpus} -i $reference -o nlrparser.txt
    """
}

process RunAssociation {
    scratch true
    cpus 1
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStats == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '6h'
    publishDir 'results', mode: 'copy'
    input:
    path presence_matrix
    path reference
    path phenotype
    path nlrparser
    output:
    path 'agrenseq_result.txt'
    script:
    """
    RunAssociation.sh -i $presence_matrix -n $nlrparser -p $phenotype -a $reference -o agrenseq_result.txt
    """
}

workflow agrenseq {
    accession_table = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) } \
        | CountKmers \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(name: "accession.tsv", newLine: true)

    phenotype_file = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> "${row.sample}\t${row.score}" } \
        | collectFile(name: "phenotype.tsv", newLine: true)

    matrix = CreatePresenceMatrix(accession_table)

    reference = Channel
        .fromPath(params.reference)

    nlrparser = NLRParser(reference)

    association = RunAssociation(matrix, reference, phenotype_file, nlrparser)
}
