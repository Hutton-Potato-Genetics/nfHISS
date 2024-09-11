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

process SamtoolsFaidx {
    //conda conda_env
    container 'swiftseal/drenseq:latest'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path reference
    output:
    path "${reference}.fai"
    script:
    """
    samtools faidx $reference
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
    cpus 1
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
        $bam
    """   
}

process IndexStrict {
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    tuple val(sample), path(strict)
    output:
    tuple val(sample), path 'strict.bam.bai'
    script:
    """
    samtools index $strict strict.bam.bai
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
        -e-value 1e-5 \
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
    container 'https://hub.docker.com/r/rocker/tidyverse/'
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
    path bait_regions_bed
    path reference_headers
    output:
    path 'passed_genes.txt'
    path 'missed_genes.txt'
    publishDir 'results', mode: 'copy'
    script:
    """
    #!/usr/bin/env python3

    bed_nlr = set()

    with open($bait_regions_bed) as bed:
        for line in bed
            bed_nlr.add(line.strip().split()[3])
    
    with open($reference_headers) as headers:
        next(headers) # skip the header
        for line in headers:
            nlr = line.strip()
            if nlr not in bed_nlr:
                with open("missing_genes.txt", "w") as missed:
                    string_to_write = nlr + " not found in bed file"
                    print(string_to_write, file = missed)
    
    with open("passed_genes.txt", "w") as passed:
        for nlr in bed_nlr:
            print(nlr, file = passed)
    """
}

process BedtoolsCoverage {
    //conda conda_env
    container 'swiftseal/drenseq:latest'
    cpus 1
    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path bed
    path bam
    path bai
    output:
    path "${bam.baseName}.coverage.txt"
    script:
    """
    bedtools coverage \
        -a $bed \
        -b $bam \
        > ${bam.baseName}.coverage.txt
    """
}

process FreeBayes {
    //conda conda_env
    container 'swiftseal/drenseq:latest'
    scratch true
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '4h'
    input:
    path reference
    path fai
    path bed
    path bam
    path bai
    output:
    tuple path("${bam.baseName}.vcf.gz"), path("${bam.baseName}.vcf.gz.tbi")
    script:
    """
    freebayes \
      -f ${reference} \
      -t ${bed} \
      --min-alternate-count 2 \
      --min-alternate-fraction 0.05 \
      --ploidy 4 \
      -m 0 \
      -v variants.vcf \
      --legacy-gls ${bam}

    bcftools sort -o variants.sorted.vcf variants.vcf

    bgzip -c variants.sorted.vcf > ${bam.baseName}.vcf.gz
    tabix -p vcf ${bam.baseName}.vcf.gz
    """
}

process MergeVCFs {
    //conda conda_env
    container 'swiftseal/drenseq:latest'
    cpus 1
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path vcf_files
    path reference
    path bed
    output:
    path 'merged.vcf'
    script:
    """
    VCF_FILES=\$(ls *.vcf.gz)
    bcftools merge -o merged.vcf \$VCF_FILES
    """
}

process CoverageMatrix{
    //container 'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '1h'
    input:
    path txt
    output:
    path 'merged.csv'
    script:
    """
    #!/usr/bin/env Rscript
    # Could use data.table, but I'm lazy :-)
    library(tidyverse)

    EXTENSION <- ".coverage.txt"

    # list all files in directory
    files <- list.files(path = "coverage", pattern = EXTENSION, full.names = TRUE)

    # read all files into a list, appending the sample name
    merged <- files %>%
      map(function(x) {
        read.table(x) %>% mutate(sample = gsub(EXTENSION, "", basename(x)))
      }) %>%
      bind_rows() %>%
      select(gene = V4, sample = sample, coverage = V8)

    matrix <- merged %>%
      pivot_wider(names_from = sample, values_from = coverage)

    write_delim(merged, "coverage_long.tsv", delim = "\t", col_names = FALSE)
    write_delim(matrix, "coverage_matrix.tsv", delim = "\t", col_names = TRUE)
    """
}

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    bed = TrimBed(file(params.bed))

    fai = SamtoolsFaidx(file(params.reference))

    reads = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) } \
        | Fastp

    (bam, bai) = BowtieAlign(bowtie2_index.first(), reads)

    (strict_bam, strict_bai) = StrictFilter(bam, bai)

    BedtoolsCoverage(bed, strict_bam, strict_bai)

    vcfs = FreeBayes(file(params.reference), fai, bed, bam, bai) \
        | collect

    MergeVCFs(vcfs, file(params.reference), bed)
}
