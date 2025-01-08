#!/usr/bin/env nextflow

process DEBREAK {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.debreak.vcf")
    
    script:
    """
    debreak --thread ${task.cpus} --min_size 30 --min_support 3 --min_quality 5 --rescue_dup --rescue_large_ins --poa --ref ${ref_fasta} --outpath ${sample}_debreak_out --bam ${bam}
    mv ${sample}_debreak_out/debreak.vcf ${sample}.${ref_fasta.simpleName}.debreak.vcf
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
    echo "${sample}" > sample.txt
    bcftools reheader -s sample.txt -o ${vcf.baseName}.bcf ${vcf}
    bcftools sort --write-index=tbi -Oz -o ${vcf}.gz ${vcf.baseName}.bcf 
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    DEBREAK(samples_ch, file(params.ref_fasta))

    BCFTOOLS_SORT_INDEX(DEBREAK.out)
}
