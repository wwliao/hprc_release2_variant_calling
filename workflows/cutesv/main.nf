#!/usr/bin/env nextflow

process CUTESV {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.cutesv.vcf")
    
    script:
    """
    cuteSV --threads ${task.cpus} --sample ${sample} --min_mapq 5 --min_support 3 --min_size 30 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 ${bam} ${ref_fasta} ${sample}.${ref_fasta.simpleName}.cutesv.vcf \${PWD}
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
    bcftools view -e 'POS==0' -Ou ${vcf} | bcftools sort --write-index=tbi -Oz -o ${vcf}.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    CUTESV(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT(CUTESV.out)
}
