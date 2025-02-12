#!/usr/bin/env nextflow

process SVDSS_SMOOTH {
    tag "${sample}"

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.smoothed.bam")

    script:
    """
    SVDSS smooth --threads ${task.cpus} --min-mapq 5 --reference ${ref_fasta} --bam ${bam} > ${sample}.${ref_fasta.baseName}.smoothed.bam
    """
}

process SAMTOOLS_INDEX {
    tag "${sample}"

    input:
    tuple val(sample), path(smoothed_bam)

    output:
    tuple val(sample), path(smoothed_bam), path("${smoothed_bam}.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${smoothed_bam} ${smoothed_bam}.bai 
    """
}

process SVDSS_SEARCH {
    tag "${sample}"

    input:
    tuple val(sample), path(smoothed_bam), path(smoothed_bai)
    path ref_fmd

    output:
    tuple val(sample), path(smoothed_bam), path(smoothed_bai), path("${sample}.${ref_fmd.baseName}.specifics.txt")

    script:
    """
    SVDSS search --threads ${task.cpus} --index ${ref_fmd} --bam ${smoothed_bam} > ${sample}.${ref_fmd.baseName}.specifics.txt
    """
}

process SVDSS_CALL {
    tag "${sample}"

    input:
    tuple val(sample), path(smoothed_bam), path(smoothed_bai), path(specifics_txt)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.svdss.vcf")

    script:
    """
    SVDSS call --threads ${task.cpus} --min-cluster-weight 3 --min-sv-length 30 --min-mapq 5 --reference ${ref_fasta} --bam ${smoothed_bam} --sfs ${specifics_txt} > ${sample}.${ref_fasta.baseName}.svdss.vcf
    """
}

process BCFTOOLS_REHEADER_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf)
    path ref_fasta

    output:
    tuple val(sample), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bcftools view --no-version -h ${vcf} > header.txt
    sed -i 's|^##reference=.*|##reference=file://${ref_fasta}|' header.txt
    sed -i 's/\tDEFAULT/\t${sample}/' header.txt
    bcftools reheader -h header.txt -o ${vcf.baseName}.reheader.vcf ${vcf}
    bcftools sort --write-index=tbi -Oz -o ${vcf}.gz ${vcf.baseName}.reheader.vcf
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    SVDSS_SMOOTH(samples_ch, file(params.ref_fasta))

    SAMTOOLS_INDEX(SVDSS_SMOOTH.out)

    SVDSS_SEARCH(SAMTOOLS_INDEX.out, file(params.ref_fmd))

    SVDSS_CALL(SVDSS_SEARCH.out, file(params.ref_fasta))

    BCFTOOLS_REHEADER_SORT(SVDSS_CALL.out, params.ref_fasta)
}
