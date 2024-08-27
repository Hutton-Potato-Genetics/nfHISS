// Trim the bed file to an expected size
process TrimBed {
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategu { task.exitStatus == 137 ? 'retry' : 'finish' }
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
    path "${sample}.bam"
    path "${sample}.bam.bai"
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
    """
}

process StrictFilter {
    //conda conda_env
    container 'swiftseal/drenseq:latest'
    scratch true
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'finish' }
    maxRetries 3
    time '2h'
    input:
    path bam
    path bai
    output:
    path "${bam.baseName}.strict.bam"
    path "${bam.baseName}.strict.bam.bai"
    script:
    """    
    sambamba view \
        --format=bam \
        --filter='[NM] == 0' \
        $bam \
        > ${bam.baseName}.strict.bam

    samtools index ${bam.baseName}.strict.bam
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
