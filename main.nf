nextflow.enable.dsl=2

/* ----------------------------------------------------
 * 1. Parameters Define
 * ---------------------------------------------------- */
params.reads      = "$projectDir/fastq_data/*_{1,2}.fastq.gz"
params.gtf        = "$projectDir/reference/Homo_sapiens.GRCh38.110.gtf"
params.hisat2_idx = "$projectDir/reference/grch38/genome"
params.outdir     = "$projectDir/results.nf"

/* ----------------------------------------------------
 * 2. Process 1: FastQC (Raw Reads QC)
 * ---------------------------------------------------- */
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    """
}

/* ----------------------------------------------------
 * 3. Process 2: fastp (Adapter & Quality Trimming)
 * ---------------------------------------------------- */
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trim_1.fastq.gz"), path("${sample_id}_trim_2.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html"
    path "${sample_id}_fastp.json"

    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_trim_1.fastq.gz -O ${sample_id}_trim_2.fastq.gz \
        -h ${sample_id}_fastp.html -j ${sample_id}_fastp.json \
        --thread 4
    """
}

/* ----------------------------------------------------
 * 4. Process 3: Alignment (HISAT2) using Trimmed Reads
 * ---------------------------------------------------- */
process HISAT2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/bams", mode: 'copy'

    input:
    tuple val(sample_id), path(trim1), path(trim2)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}_sorted.bam.bai"), emit: bai

    script:
    """
    hisat2 -p 2 -x ${params.hisat2_idx} -1 ${trim1} -2 ${trim2} | \
    samtools view -bS - | \
    samtools sort -o ${sample_id}_sorted.bam -
    samtools index ${sample_id}_sorted.bam
    """
}

/* ----------------------------------------------------
 * 5. Process 4: FeatureCounts (Matrix Generation)
 * ---------------------------------------------------- */
process FEATURECOUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path bams

    output:
    path "gene_counts.txt"

    script:
    """
    featureCounts -p -t exon -g gene_id \
        -a ${params.gtf} \
        -o gene_counts.txt \
        -T 4 \
        ${bams}
    """
}

/* ----------------------------------------------------
 * 6. Master Workflow Execution
 * ---------------------------------------------------- */
workflow {
    // 1. Read input FASTQ pairs
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // 2. Raw QC
    FASTQC(read_pairs_ch)

    // 3. Trimming
    FASTP(read_pairs_ch)

    // 4. Alignment using TRIMMED reads output from FASTP
    HISAT2_ALIGN(FASTP.out.trimmed_reads)

    // 5. Collect BAMs & Count
    all_bams_ch = HISAT2_ALIGN.out.bam.map{ sample_id, bam -> bam }.collect()
    FEATURECOUNTS(all_bams_ch)
}

