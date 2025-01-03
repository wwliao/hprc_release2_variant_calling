#!/usr/bin/env nextflow

process SVIM_ALIGNMENT {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.svim.vcf")
    
    script:
    """
    svim alignment --min_mapq 5 --min_sv_size 30 --minimum_depth 3 --sample ${sample} --interspersed_duplications_as_insertions ${sample}_svim_results ${bam} ${ref_fasta}
    mv ${sample}_svim_results/variants.vcf ${sample}.${ref_fasta.simpleName}.svim.vcf
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
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    SVIM_ALIGNMENT(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT_INDEX(SVIM_ALIGNMENT.out)
}
