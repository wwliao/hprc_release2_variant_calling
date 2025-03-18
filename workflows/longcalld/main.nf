#!/usr/bin/env nextflow

process LONGCALLD_CALL {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.longcalld.vcf"), emit: vcf
    tuple val(sample), path("${bam.baseName.tokenize('.')[0..1].join('.')}.longcalld.phased.bam"), emit: bam

    script:
    """
    longcallD call --threads ${task.cpus} --hifi --min-cov 6 --alt-cov 3 --min-mapq 5 --sample-name ${sample} --out-bam ${bam.baseName.tokenize('.')[0..1].join('.')}.longcalld.phased.bam --out-vcf ${sample}.${ref_fasta.baseName}.longcalld.vcf ${ref_fasta} ${bam}
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

process SAMTOOLS_EXTRACT {
    tag "${sample}"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${bam.baseName}.tsv.gz")

    script:
    """
    (
        echo -e "READ_NAME\tCHROM\tPOS\tHP\tPS"
        samtools view ${bam} | awk '
        {
            hp = "."; ps = ".";
            for(i=12; i<=NF; i++) {
                if(\$i ~ /^HP:/) hp = substr(\$i, 6);
                if(\$i ~ /^PS:/) ps = substr(\$i, 6);
            }
            print \$1 "\t" \$3 "\t" \$4 "\t" hp "\t" ps;
        }'
    ) | gzip > ${bam.baseName}.tsv.gz
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam)) }
        .set { samples_ch }
    
    LONGCALLD_CALL(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_SORT(LONGCALLD_CALL.out.vcf)

    SAMTOOLS_EXTRACT(LONGCALLD_CALL.out.bam)
}
