conda_env = 'bioconda::fastp=0.23.4 bioconda::freebayes=1.3.7 bioconda::sambamba=1.0 bioconda::samtools=1.19.2 bioconda::bcftools=1.19 bioconda::bowtie2=2.5.3 bioconda::htslib=1.19.1 bioconda::bedtools=2.31.1'

process Fastp {
    conda conda_env
    scratch true
    cpus 1
    memory { 4.GB * task.attempt }
    errorStrategy {'retry'}
    maxRetries 3
    time '4h'
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path('R1.fastq.gz'), path('R2.fastq.gz')
    script:
    """
    micromamba run fastp -i $read1 -I $read2 -o R1.fastq.gz -O R2.fastq.gz
    """
}


process BowtieBuild {
    conda conda_env
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
    micromamba run bowtie2-build $reference bowtie2_index/reference
    """
}

process BowtieAlign {
    conda conda_env
    scratch true
    cpus 8
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    time '4h'
    input:
    path bowtie2_index
    tuple val(sample), path(read1), path(read2)
    output:
    path "${sample}.bam"
    path "${sample}.bam.bai"
    script:
    """
    micromamba run bowtie2 \
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
    
    micromamba run sambamba view \
        --format=bam \
        --filter='[NM] == 0' \
        aligned.bam \
        > ${sample}.bam

    micromamba run samtools index ${sample}.bam
    """
}

process BedtoolsCoverage {
    publishDir 'coverage', mode: 'copy'
    conda conda_env
    cpus 1
    memory { 1.GB * task.attempt }
    errorStrategy 'retry'
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
    micromamba run bedtools coverage \
        -a $bed \
        -b $bam \
        > ${bam.baseName}.coverage.txt
    """
}

process FreeBayes {
    conda conda_env
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
    path bai
    output:
    tuple path("${bam.baseName}.vcf.gz"), path("${bam.baseName}.vcf.gz.tbi")
    script:
    """
    micromamba run freebayes \
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

    micromamba run bcftools sort -o variants.sorted.vcf variants.vcf

    micromamba run bgzip -c variants.sorted.vcf > ${bam.baseName}.vcf.gz
    micromamba run tabix -p vcf ${bam.baseName}.vcf.gz
    """
}

process MergeVCFs {
    conda conda_env
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
    VCF_FILES=\$(ls *.vcf.gz)
    micromamba run bcftools merge -o merged.vcf $VCF_FILES
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

    (bam, bai) = BowtieAlign(bowtie2_index.first(), reads)

    BedtoolsCoverage(bed, bam, bai)

    vcfs = FreeBayes(file(params.reference), bed, bam, bai) \
        | collect

    MergeVCFs(vcfs, file(params.reference), bed)
}
