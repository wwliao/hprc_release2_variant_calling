#!/usr/bin/env nextflow

process KANPIG_GT {
    tag "${sample}"

    input:
    tuple val(sample), val(sex), path(query_bam), path(query_bai), path(query_vcf)
    path ref_fasta
    path ref_fai
    path ploidy_male_bed
    path ploidy_female_bed

    output:
    tuple val(sample), path("${query_vcf.baseName.tokenize('.')[0..2].join('.')}.genotyped.vcf")

    script:
    if( sex == 'female' )
        """
        kanpig gt --threads ${task.cpus} --input ${query_vcf} --reads ${query_bam} --reference ${ref_fasta} --ploidy-bed ${ploidy_female_bed} --sizemin 30 --out ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.genotyped.vcf
        """

    else if( sex == 'male' )
        """
        kanpig gt --threads ${task.cpus} --input ${query_vcf} --reads ${query_bam} --reference ${ref_fasta} --ploidy-bed ${ploidy_male_bed} --sizemin 30 --out ${query_vcf.baseName.tokenize('.')[0..2].join('.')}.genotyped.vcf
        """
}

process BCFTOOLS_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(genotyped_vcf)

    output:
    tuple val(sample), path("${genotyped_vcf.baseName}.vcf.gz"), path("${genotyped_vcf.baseName}.vcf.gz.tbi")

    script:
    """
    bcftools sort --write-index=tbi -Oz -o ${genotyped_vcf.baseName}.vcf.gz ${genotyped_vcf}
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.sex, file(row.query_bam), file(row.query_bai), file(row.query_vcf)) }
        .set { samples_ch }

    KANPIG_GT(samples_ch, file(params.ref_fasta), file(params.ref_fai), file(params.ploidy_male_bed), file(params.ploidy_female_bed))

    BCFTOOLS_SORT(KANPIG_GT.out)
}
