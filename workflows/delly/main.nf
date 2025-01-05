#!/usr/bin/env nextflow

process DELLY_LR {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.delly.bcf")
    
    script:
    """
    delly lr --technology pb --mapqual 5 --min-clique-size 3 --geno-qual 5 --genome ${ref_fasta} --outfile ${sample}.${ref_fasta.simpleName}.delly.bcf ${bam}
    """
}

process BCFTOOLS_SORT_INDEX {
    tag "${sample}"

    input:
    tuple val(sample), path(bcf)

    output:
    tuple val(sample), path("${bcf.baseName}.vcf.gz"), path("${bcf.baseName}.vcf.gz.tbi")

    script:
    """
    bcftools sort --write-index=tbi -Oz -o ${bcf.baseName}.vcf.gz ${bcf}
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    DELLY_LR(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT_INDEX(DELLY_LR.out)
}
