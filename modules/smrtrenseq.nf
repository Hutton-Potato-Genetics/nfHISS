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

process CanuAssemble {
    container 'docker://quay.io/biocontainers/canu:2.2--ha47f30e_0'
    scratch true
    cpus 8
    memory { 36.GB * task.attempt }
    maxRetries 3
    time '48h'
    input:
    tuple val(sample), path(reads)
    val genome_size
    val max_input_coverage
    output:
    path "assembly/${sample}_assembly.contigs.fasta"
    path "assembly/${sample}.report"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    canu -d assembly -p assembly genomeSize=$genome_size useGrid=false -pacbio-hifi $reads maxInputCoverage=$max_input_coverage batMemory=32g
    """
}

process SeqkitStats {
    container 'docker://quay.io/biocontainers/seqfu:1.20.3--h1eb128b_2'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    maxRetries 3
    time '1h'
    input:
    path assembly
    tuple val(sample), path(reads)
    output:
    path "${sample}_statistics.txt"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    seqkit stats -b $assembly | sed 's/_assembly\.contigs//g' > ${sample}_statistics.txt
    """
}

process ChopSequences {
    container 'docker://quay.io/biocontainers/meme:5.5.6--pl5321h4242488_0'
    scratch true
    cpus 1
    memory { 2.GB * task.attempt }
    maxRetries 3
    time '2h'
    input:
    path assembly
    output:
    path "chopped.fa"
    script:
    """
    chop_sequences.sh -i $assembly -o chopped.fa
    """
}

process NLRParser {
    container 'docker://quay.io/biocontainers/meme:5.5.6--pl5321h4242488_0'
    scratch true
    cpus 2
    memory { 4.GB * task.attempt }
    maxRetries 3
    time '8h'
    input:
    path chopped
    output:
    path "parser.xml"
    script:
    """
    nlr_parser.sh -t 2 -i $chopped -o parser.xml
    """
}

process NLRAnnotator {
    container 'docker://quay.io/biocontainers/meme:5.5.6--pl5321h4242488_0'
    scratch true
    cpus 1
    memory { 2.GB * task.attempt }
    maxRetries 3
    time '4h'
    input:
    path assembly
    path parser_xml
    tuple val(sample), path(reads)
    val flanking
    output:
    path "${sample}_NLR_annotator.txt"
    path "${sample}_NLR_annotator.fa"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    nlr_annotator.sh -i $parser_xml -o ${sample}_NLR_annotator.txt -f $assembly ${sample}_NLR_annotator.fa $flanking
    """
}

process SummariseNLRs {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    maxRetries 3
    time '2h'
    input:
    path annotator_text
    tuple val(sample), path(reads)
    output:
    path "${sample}_NLR_summary.txt"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    summarise_nlrs.py --input $annotator_text --output ${sample}_NLR_summary.txt
    """
}

process InputStatistics {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    maxRetries 3
    time '2h'
    input:
    path report
    tuple val(sample), path(reads)
    output:
    path "${sample}_input_stats.txt"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    input_stats.sh $report ${sample}_input_stats.txt
    """
}

process NLR2Bed {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    maxRetries 3
    time '2h'
    input:
    path annotator_text
    tuple val(sample), path(reads)
    output:
    path "${sample}_NLR_Annotator.bed"
    publishDir "results/${sample}", mode: 'copy'
    script:
    """
    nlr2bed.py --input $annotator_text --output ${sample}_NLR_Annotator.bed
    """
}

workflow smrtrenseq {
    trimmed_reads = Cutadapt(reads)
    assembly = Canu(trimmed_reads)
}
