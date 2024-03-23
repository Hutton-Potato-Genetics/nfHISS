
process Fastp {
    conda 'envs/drenseq.yml'
    scratch true
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path('R1.fastq.gz'), path('R2.fastq.gz')
    script:
    """
    fastp -i $read1 -I $read2 -o R1.fastq.gz -O R2.fastq.gz
    """
}


process BowtieBuild {
    conda 'envs/drenseq.yml'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
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
    conda 'envs/drenseq.yml'
    scratch true
    cpus 8
    memory { 2.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    path bowtie2_index
    tuple val(sample), path(read1), path(read2)
    output:
    path "${sample}.bam"
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
      | samtools sort -@ ${task.cpus} -o aligned.bam
    
    sambamba view \
        --format=bam \
        --filter='[NM] == 0' \
        aligned.bam \
        > ${sample}.bam
    """
}

process BedtoolsCoverage {
    publishDir 'coverage', mode: 'copy'
    conda 'envs/drenseq.yml'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '1h'
    input:
    path bed
    path bam
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
    conda 'envs/drenseq.yml'
    scratch true
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    path reference
    path bed
    path bam
    output:
    tuple path("${vcf.baseName}.vcf.gz"), path("${vcf.baseName}.vcf.gz.tbi")
    script:
    """
    freebayes \
      -f ${reference} \
      -t ${bed} \
      --min-alternate-count 2 \
      --min-alternate-fraction 0.05 \
      --ploidy 2 \
      --throw-away-indel-obs \
      --throw-away-mnps-obs \
      --throw-away-complex-obs \
      -m 0 \
      -v variants.vcf \
      --legacy-gls ${bam}

    bcftools sort -o variants.sorted.vcf variants.vcf

    bgzip -c variants.sorted.vcf > ${bam.baseName}.vcf.gz
    tabix -p vcf ${bam.baseName}.vcf.gz
    """
}

process MergeVCFs {
    conda 'envs/drenseq.yml'
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
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
    VCF_FILES=\$()
    bcftools merge -o merged.vcf $VCF_FILES
    """
}

workflow drenseq {
    bowtie2_index = Channel
        .fromPath(params.reference) \
        | BowtieBuild

    bed = file(params.bed)

    reads = Channel
        .fromPath(params.reads)
        .splitCsv(header: true, sep: "\t")
        .map { row -> tuple(row.sample, file(row.forward), file(row.reverse)) } \
        | Fastp

    bams = BowtieAlign(bowtie2_index.first(), reads)

    BedtoolsCoverage(bed, bams)

    FreeBayes(file(params.reference), bed, bams) \
        | collect \
        | MergeVCFs
}
