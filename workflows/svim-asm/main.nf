#!/usr/bin/env nextflow

process SVIMASM_DIPLOID {
    tag "${sample}"

    input:
    tuple val(sample), path(bam1), path(bai1), path(bam2), path(bai2)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.svim-asm.vcf")

    script:
    """
    svim-asm diploid --min_mapq 5 --min_sv_size 30 --sample ${sample} --query_names --interspersed_duplications_as_insertions ${sample}_svimasm_results ${bam1} ${bam2} ${ref_fasta}
    mv ${sample}_svimasm_results/variants.vcf ${sample}.${ref_fasta.baseName}.svim-asm.vcf
    """
}

process BCFTOOLS_SORT_INDEX {
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
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam1), file(row.bai1), file(row.bam2), file(row.bai2)) }
        .set { samples_ch }

    SVIMASM_DIPLOID(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT_INDEX(SVIMASM_DIPLOID.out)
}
