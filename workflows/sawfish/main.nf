#!/usr/bin/env nextflow

process SAWFISH {
    tag "${sample}"
    
    input:
    tuple val(sample), val(sex), path(bam), path(bai)
    path ref_fasta
    path ref_cn_xx
    path ref_cn_xy
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz"), path("${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz.tbi")
    
    script:
    if( sex == 'female' )
        """
        sawfish discover --threads ${task.cpus} --min-indel-size 30 --min-sv-mapq 5 --disable-path-canonicalization --ref ${ref_fasta} --expected-cn ${ref_cn_xx} --bam ${bam} --output-dir ${sample}_sawfish_discover_output && \
        sawfish joint-call --threads ${task.cpus} --min-sv-mapq 5 --sample ${sample}_sawfish_discover_output --output-dir ${sample}_sawfish_joint-call_output
        mv ${sample}_sawfish_joint-call_output/genotyped.sv.vcf.gz ${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz
        mv ${sample}_sawfish_joint-call_output/genotyped.sv.vcf.gz.tbi ${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz.tbi
        """

    else if( sex == 'male' )
        """
        sawfish discover --threads ${task.cpus} --min-indel-size 30 --min-sv-mapq 5 --disable-path-canonicalization --ref ${ref_fasta} --expected-cn ${ref_cn_xy} --bam ${bam} --output-dir ${sample}_sawfish_discover_output && \
        sawfish joint-call --threads ${task.cpus} --min-sv-mapq 5 --sample ${sample}_sawfish_discover_output --output-dir ${sample}_sawfish_joint-call_output
        mv ${sample}_sawfish_joint-call_output/genotyped.sv.vcf.gz ${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz
        mv ${sample}_sawfish_joint-call_output/genotyped.sv.vcf.gz.tbi ${sample}.${ref_fasta.simpleName}.sawfish.vcf.gz.tbi
        """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.sex, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    SAWFISH(samples_ch, file(params.ref_fasta), file(params.ref_cn_xx), file(params.ref_cn_xy))
}
