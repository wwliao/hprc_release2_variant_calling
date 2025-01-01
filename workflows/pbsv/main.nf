#!/usr/bin/env nextflow

process PBSV_DISCOVER {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam)
    path ref_trf
    
    output:
    tuple val(sample), path("${sample}.${ref_trf.simpleName}.svsig.gz")
    
    script:
    """
    pbsv discover --hifi --sample ${sample} --tandem-repeats ${ref_trf} ${bam} ${sample}.${ref_trf.simpleName}.svsig.gz
    """
}

process PBSV_CALL {
    tag "${sample}"
    
    input:
    tuple val(sample), path(svsig)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.pbsv.vcf")
    
    script:
    """
    pbsv call --hifi --preserve-non-acgt --min-sv-length 40 --num-threads ${task.cpus} ${ref_fasta} ${svsig} ${sample}.${ref_fasta.simpleName}.pbsv.vcf
    """
}

process BCFTOOLS_SORT_INDEX {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.pbsv.vcf.gz"), path("${sample}.${ref_fasta.simpleName}.pbsv.vcf.gz.tbi")

    script:
    """
    bcftools sort --write-index=tbi -Oz -o ${sample}.${ref_fasta.simpleName}.pbsv.vcf.gz ${vcf} 
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam)) }
        .set { samples_ch }

    PBSV_DISCOVER(samples_ch, file(params.ref_trf))

    PBSV_CALL(PBSV_DISCOVER.out, file(params.ref_fasta))

    BCFTOOLS_SORT_INDEX(PBSV_CALL.out, file(params.ref_fasta))
}
