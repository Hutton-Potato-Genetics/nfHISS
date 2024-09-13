process TrimBed {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    time '1h'
    input:
    path bed
    output:
    path 'trimmed.bed'
    script:
    """
    awk '{print \$1"\t"\$2"\t"\$3"\t"\$4}' $bed > trimmed.bed
    """
}

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

process BowtieBuild {
    container 'docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
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
    container 'docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4'
    scratch true
    cpus 8
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '6h'
    input:
    path bowtie2_index
    tuple val(sample), path(read1), path(read2)
    val score
    val max_align
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
      --score-min $score \
      --phred33 \
      --fr \
      --maxins 1000 \
      --very-sensitive \
      --no-unal \
      --no-discordant \
      -k $max_align \
      -S aligned.sam
    """
}

process ParseAlignment {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'
    scratch true
    cpus 2
    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path('aligned.bam')
    tuple val(sample), path('aligned.bam.bai')
    script:
    """
    samtools view $sam -b -o aligned_unsorted.bam -@ $task.cpus
    samtools sort -l 9 $sam -o aligned.bam -@ $task.cpus
    samtools index aligned.bam aligned.bam.bai -@ $task.cpus
    """
}

process StrictFilter {
    container 'docker://quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_2'
    scratch true
    cpus 2
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    tuple val(sample), path(bam)
    tuple val(sample), path(bai)
    output:
    tuple val(sample), path('strict.bam')
    script:
    """    
    sambamba view \
        --format=bam \
        -l 9 \
        --filter='[NM] == 0' \
        -o strict.bam \
        -t 2 \
        $bam
    """   
}

process BaitsBlasting {
    container 'docker://quay.io/biocontainers/blast:2.16.0--hc155240_2'
    scratch true
    cpus 8
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '8h'
    input:
    path reference
    path baits
    val identity
    val coverage
    output:
    path 'blast_out.txt'
    script:
    """
    makeblastdb -in $reference -out blast_db -dbtype nucl
    blastn \
        -db blast_db \
        -query $baits \
        -out blast_out.txt \
        -perc_identity $identity \
        -qcov_hsp_perc $coverage \
        -evalue 1e-5 \
        -num_threads $task.cpus
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen len qcovs qcovhsp'
    """
}

process Headers {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path trimmed_bed
    path reference
    output:
    path 'nlr_headers.txt'
    path 'reference_headers.txt'
    script:
    """
    echo "gene" > nlr_headers.txt
    cat $trimmed_bed | cut -f4 >> nlr_headers.txt
    cat $reference | grep '>' | sed 's/>//g' > reference_headers.txt
    """
}

process IdentifyBaitRegions {
    container 'docker://quay.io/biocontainers/bioconductor-biostrings:2.70.1--r43ha9d7317_2'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    path blast_out
    path reference_headers
    path reference
    val flank
    output:
    path 'bait_regions.bed'
    script:
    """
    RangeReduction.R $blast_out bait_regions.bed $flank $reference_headers $reference
    """
}

process AnnotatorBaits {
    container 'docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    path bait_regions
    path trimmed_bed
    output:
    path 'nlr_bait_regions.bed'
    script:
    """
    bedtools intersect -a $trimmed_bed -b $bait_regions > nlr_bait_regions.bed
    """
}

process BaitBlastCheck {
    container 'docker://quay.io/biocontainers/python:3.12'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path nlr_bait_regions_bed
    path reference_headers
    output:
    path 'passed_genes.txt'
    path 'missing_genes.txt'
    publishDir 'results/diagnostics', mode: 'copy'
    script:
    """
    #!/usr/bin/env python3

    bed_nlr = set()

    with open("$nlr_bait_regions_bed") as bed:
        for line in bed:
            bed_nlr.add(line.strip().split()[3])
    
    with open("$reference_headers") as headers:
        next(headers) # skip the header
        for line in headers:
            with open("missing_genes.txt", "w") as missed:
                nlr = line.strip()
                if nlr not in bed_nlr:
                    string_to_write = nlr + " not found in bed file"
                    print(string_to_write, file = missed)
        with open("missing_genes.txt", "w") as missed:
            print("\\n", file = missed)
    
    with open("passed_genes.txt", "w") as passed:
        for nlr in bed_nlr:
            print(nlr, file = passed)
    """
}

process BedtoolsCoverage {
    container 'docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path bait_regions_bed
    tuple val(sample), path(bam)
    path passed_genes
    output:
    tuple val(sample), path('coverage.txt')
    script:
    """
    bedtools coverage -d \
        -a $bait_regions_bed \
        -b $bam \
        > coverage.txt
    """
}

process PerGeneCoverage {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path nlr_headers
    tuple val(sample), path(coverage)
    output:
    tuple val(sample), path('gene_coverage.txt')
    script:
    """
    cat $nlr_headers | tail -n +2 | while read gene
    do
        numPosWithCoverage=`grep -w "\$gene" $coverage | awk '\$6>0' | wc -l`
        numPosTotal=`grep -w "\$gene" $coverage | wc -l` 
        if [ \$numPosTotal -eq 0 ]
        then
            echo "ERROR: gene \$gene has CDS region of length zero. Check your input data (e.g. gene spelling in FASTA and CDS BED file) and retry.\nAborting pipeline run."
            exit 1
        fi
        pctCov=`awk "BEGIN {{print (\$numPosWithCoverage/\$numPosTotal)*100 }}"`
        echo -e "\n# covered positions for sample $sample in gene \$gene: \$numPosWithCoverage\n# CDS positions for gene \$gene: \$numPosTotal\npctCov: \$pctCov"
        echo -e "\$gene\t\$pctCov" >> gene_coverage.txt
    done
    """
}

process CombineGeneCoverages {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    tuple val(sample), path(gene_coverage)
    output:
    path "${sample}_coverage_values.txt"
    script:
    """
    echo $sample > ${sample}_coverage_values.txt
    cat $gene_coverage | cut -f2 >> ${sample}_coverage_values.txt
    """
}

process CombineCoverageValues {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path coverage_values
    path nlr_headers
    val ulimit
    output:
    path 'all_coverage_values.txt'
    script:
    """
    ulimit -n $ulimit
    paste $nlr_headers $coverage_values > all_coverage_values.txt
    """
}

process TransposeCombinedCoverage {
    container 'docker://quay.io/biocontainers/pandas:2.2.1'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path all_coverage_values
    output:
    path 'all_coverage_values_transposed.txt'
    publishDir 'results', mode: 'copy'
    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd

    df = pd.read_table("$all_coverage_values", header = None)
    df.T.to_csv("all_coverage_values_transposed.txt", sep = "\t", header = False, index = False)
    """
}

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    trimmed_bed = TrimBed(file(params.bed))

    reads = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.FRead), file(row.RRead)) }
    
    trimmed_reads = TrimReads(reads, params.adaptor_1, params.adaptor_2)

    sam = BowtieAlign(bowtie2_index.first(), trimmed_reads, params.score, params.max_align)

    (bam, bai) = ParseAlignment(sam)

    strict_bam = StrictFilter(bam, bai)

    blast_out = BaitsBlasting(params.reference, params.baits, params.identity, params.coverage)

    (nlr_headers, reference_headers) = Headers(trimmed_bed, params.reference)

    bait_regions_bed = IdentifyBaitRegions(blast_out, reference_headers, params.reference, params.flank)

    nlr_bait_regions_bed = AnnotatorBaits(bait_regions_bed, trimmed_bed)

    (passed, missed) = BaitBlastCheck(nlr_bait_regions_bed, reference_headers)

    coverage = BedtoolsCoverage(nlr_bait_regions_bed, strict_bam, passed)

    gene_coverage = PerGeneCoverage(nlr_headers, coverage)

    sample_coverage = CombineGeneCoverages(gene_coverage)

    all_coverage_values = CombineCoverageValues(sample_coverage.collect(), nlr_headers, params.ulimit)

    transposed_coverage = TransposeCombinedCoverage(all_coverage_values)
}
