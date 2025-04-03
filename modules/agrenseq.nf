process TrimReads {
    container 'docker://quay.io/biocontainers/cutadapt:4.9--py312hf67a6ed_0'
    scratch true
    cpus 8
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
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
        -a file:$adaptor_1 \
        -A file:$adaptor_2 \
        -o R1.fastq.gz \
        -p R2.fastq.gz \
        $read1 \
        $read2
    """
}

process CountKmers {
    container 'docker://quay.io/biocontainers/kmc:3.2.4--h6dccd9a_1'
    scratch true
    cpus 2
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path("${sample}.dump")
    script:
    """
    cat $read1 $read2 > reads.fq.gz
    kmc -k51 -m${task.memory.toGiga()} -t${task.cpus} reads.fq.gz kmc_output .
    kmc_tools -t${task.cpus} transform kmc_output -ci10 dump ${sample}.dump
    """
}

process CreatePresenceMatrix {
    scratch true
    cpus 1
    memory { 32.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 1.h * task.attempt }
    input:
    path accession_table
    output:
    path 'presence_matrix.txt'
    script:
    """
    create_presence_matrix.sh ${task.memory.toGiga()}G -i $accession_table -o presence_matrix.txt
    """
}

process NLRParser {
    container 'community.wave.seqera.io/library/meme_openjdk:3e840cb4617be872'
    scratch true
    cpus 4
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
    input:
    path reference
    output:
    path 'nlrparser.txt'
    script:
    """
    nlr_parser.sh -Xmx${task.memory.toMega()}M -t ${task.cpus} -i $reference -o nlrparser.txt
    """
}

process RunAssociation {
    scratch true
    cpus 2
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 15.m * task.attempt }
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
    run_association.sh ${task.memory.toGiga()}G -i $presence_matrix -n $nlrparser -p $phenotype -a $reference -o agrenseq_result.txt
    """
}

process Blast {
    container 'docker://quay.io/biocontainers/blast:2.16.0--hc155240_2'
    scratch true
    cpus 8
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 45.m * task.attempt }
    input:
    path blast_reference
    path association_reference
    output:
    path 'blast.txt'
    script:
    """
    makeblastdb -in $blast_reference -dbtype nucl -out blast_ref
    blastn -query $association_reference -db blast_ref -outfmt 6 -num_threads $task.cpus > blast.txt
    """
}

process SortBlast {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
    input:
    path blast_out
    output:
    path 'blast_sorted.txt'
    script:
    """
    cat $blast_out | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 -m > blast_sorted.txt
    """
}

process GetSizes {
    container 'docker://quay.io/biocontainers/bioawk:1.0--he4a0461_12'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
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
    container 'docker://rocker/tidyverse:4.4.1'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 10.m * task.attempt }
    input:
    path blast_text
    path sizes
    path association_results
    val threshold
    val title
    output:
    path 'AgRenSeq_plot.png'
    path 'Blast_plot.png'
    path 'filtered_contigs.txt'
    publishDir 'results', mode: 'copy'
    script:
    """
    plot.R $association_results $threshold $title filtered_contigs.txt AgRenSeq_plot.png
    blast_plot.R $sizes $blast_text filtered_contigs.txt $title Blast_plot.png
    """
}

process FinalFilePrep {
    scratch true
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    maxRetries 3
    time { 15.m * task.attempt }
    input:
    path association_reference
    path annotator_bed
    path filtered_contigs
    output:
    path 'candidates.fa'
    path 'candidates.bed'
    publishDir 'results', mode: 'copy'
    script:
    """
    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < $association_reference | tail -n +2 > unwrapped.fa
    cat unwrapped.fa | grep -A1 -f $filtered_contigs | sed 's/--//g' | sed '/^\$/d' | sed '/^>/ s/ .*//' > candidates.fa
    cat $annotator_bed | grep -f $filtered_contigs | cut -f1-4 > candidates.bed
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

    blast_raw = Blast(params.blast_reference, association_reference)

    blast_text = SortBlast(blast_raw)

    sizes = GetSizes(params.blast_reference)

    (ag_plot, blast_plot, filtered_contigs) = Plot(blast_text, sizes, association, params.threshold, params.title)

    (candidates_fa, candidates_bed) = FinalFilePrep(association_reference, params.annotator_bed, filtered_contigs)
}
