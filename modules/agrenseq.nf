process TrimReads {
    container 'docker://quay.io/biocontainers/cutadapt:4.9--py312hf67a6ed_0'
    scratch true
    cpus 8
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(read1), path(read2)
    path adaptor_1
    path adaptor_2
    output:
    tuple val(sample), path('R1.fastq.gz'), path('R2.fastq.gz')
    script:
    """
    cutadapt \
        -j $task.cpus \
        --minimum-length 50 \
        -q 20,20 \
        -a file:$adaptor_1 \
        -A file:$adaptor_2 \
        -o R1.fastq.gz \
        -p R2.fastq.gz \
        $read1 \
        $read2
    """
}

process CountKmers {
    container 'docker://quay.io/biocontainers/kmc3.2.4--h6dccd9a_1'
    scratch true
    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
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
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    path accession_table
    output:
    path 'presence_matrix.txt'
    script:
    """
    create_presence_matrix.sh -i $accession_table -o presence_matrix.txt
    """
}

process NLRParser {
    conda 'bioconda::meme=5.4.1 conda-forge::openjdk=11.0.1'
    scratch true
    cpus 4
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
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
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '6h'
    input:
    path presence_matrix
    path reference
    path phenotype
    path nlrparser
    output:
    path 'agrenseq_result.txt'
    publishDir 'results', mode: 'copy'
    script:
    """
    run_association.sh -i $presence_matrix -n $nlrparser -p $phenotype -a $reference -o agrenseq_result.txt
    """
}

process Blast {
    container 'docker://quay.io/biocontainers/blast:2.16.0--hc155240_2'
    scratch true
    cpus 8
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '8h'
    input:
    path blast_reference
    path association_reference
    output:
    path 'blast_sorted.txt'
    script:
    """
    makeblastdb -in $blast_reference -dbtype nucl -out blast_ref
    blastn -query $association_reference -db blast_ref -outfmt 6 -num_threads $task.cpus | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > blast_sorted.txt
    """
}

process GetSizes {
    container 'docker://quay.io/biocontainers/bioawk:1.0--he4a0461_12'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path blast_reference
    output:
    path 'sizes.txt'
    script:
    """
    bioawk -c fastx '{{ print \$name, length(\$seq) }}' $blast_reference > sizes.txt
    """
}

process Plot {
    container 'https://hub.docker.com/r/rocker/tidyverse/'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '4h'
    input:
    path blast_text
    path sizes
    path association_results
    val threshold
    val title
    output:
    path 'AgRenSeq_plot.png'
    path 'Blast_plot.png'
    script:
    """
    plot.R $association_results $threshold $title filtered_contigs.txt AgRenSeq_plot.png
    blast_plot.R $sizes $blast_text filtered_contigs.txt $title Blast_plot.png
    """
}

workflow agrenseq {
    reads = Channel.fromPath(params.reads).splitCsv(header: true, sep: "\t").map { row -> tuple(row.sample, file(row.R1), file(row.R2)) }

    trimmed_reads = TrimReads(reads, params.adaptor_1, params.adaptor_2)

    accession_table =  CountKmers(trimmed_reads) \
        | map { it -> "${it[0]}\t${it[1]}" } \
        | collectFile(name: "accession.tsv", newLine: true)

    phenotype_file = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> "${row.sample}\t${row.score}" } \
        | collectFile(name: "phenotype.tsv", newLine: true)

    matrix = CreatePresenceMatrix(accession_table)

    association_reference = Channel
        .fromPath(params.association_reference)

    nlrparser = NLRParser(association_reference)

    association = RunAssociation(matrix, association_reference, phenotype_file, nlrparser)

    blast_text = Blast(params.blast_reference, association_reference)

    sizes = GetSizes(params.blast_reference)

    (ag_plot, blast_plot) = Plot(blast_text, sizes, association, params.threshold, params.title)
}
