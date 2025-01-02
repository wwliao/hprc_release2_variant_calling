#!/usr/bin/env nextflow

process SNIFFLES {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_trf
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.sniffles.vcf.gz"), path("${sample}.${ref_fasta.simpleName}.sniffles.vcf.gz.tbi"), path("${sample}.${ref_fasta.simpleName}.snf")
    
    script:
    """
    sniffles --minsupport 3 --minsvlen 30 --mapq 5 --threads ${task.cpus} --reference ${ref_fasta} --tandem-repeats ${ref_trf} --sample-id ${sample} --input ${bam} --vcf ${sample}.${ref_fasta.simpleName}.sniffles.vcf.gz --snf ${sample}.${ref_fasta.simpleName}.snf
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    SNIFFLES(samples_ch, file(params.ref_fasta), file(params.ref_trf))
}
