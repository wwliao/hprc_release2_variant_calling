#!/usr/bin/env nextflow

process SAMTOOLS_MERGE_SORT {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam1), path(bam2)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.sorted.bam"), path("${sample}.${ref_fasta.baseName}.sorted.bam.bai")
    
    script:
    """
    samtools merge -@ ${task.cpus} - ${bam1} ${bam2} | samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sample}.${ref_fasta.baseName}.sorted.bam##idx##${sample}.${ref_fasta.baseName}.sorted.bam.bai -
    """
}

process CUTESVASM {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.cutesv-asm.vcf")
    
    script:
    """
    cuteSV --threads ${task.cpus} --sample ${sample} --min_mapq 5 --min_support 1 --min_size 30 --max_size -1 --report_readid --genotype --max_split_parts -1 --merge_del_threshold 500 --merge_ins_threshold 500 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 ${bam} ${ref_fasta} ${sample}.${ref_fasta.simpleName}.cutesv-asm.intermediate.vcf \${PWD} && \
    /opt/cutesv-asm/src/cuteSV/diploid_calling.py --hap1 '${sample}#1' --hap2 '${sample}#2' ${sample}.${ref_fasta.simpleName}.cutesv-asm.intermediate.vcf ${sample}.${ref_fasta.simpleName}.cutesv-asm.vcf
    """
}

process BCFTOOLS_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf)

    output:
    tuple val(sample), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bcftools sort --write-index=tbi -Oz -o ${vcf}.gz ${vcf} 
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam1), file(row.bam2)) }
        .set { samples_ch }

    SAMTOOLS_MERGE_SORT(samples_ch, file(params.ref_fasta))

    CUTESVASM(SAMTOOLS_MERGE_SORT.out, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT(CUTESVASM.out)
}
