#!/usr/bin/env nextflow

process SVISIONPRO {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.svision_pro_v2.4.s3.vcf")

    script:
    """
    SVision-pro --process_num ${task.cpus} --preset hifi --min_supp 3 --min_mapq 5 --min_sv_size 30 --rescue_large --detect_mode germline --img_size 512 --model_path /opt/svision-pro/src/pre_process/model_liteunet_256_8_16_32_32_32.pth --sample_name ${sample} --target_path ${bam} --genome_path ${ref_fasta} --out_path \${PWD}
    """
}

process BCFTOOLS_REHEADER_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.svision-pro.vcf.gz"), path("${sample}.${ref_fasta.baseName}.svision-pro.vcf.gz.tbi")

    script:
    """
    echo "${sample}" > sample.txt
    bcftools reheader -s sample.txt -o ${sample}.${ref_fasta.baseName}.svision-pro.bcf ${vcf}
    bcftools sort --write-index=tbi -Oz -o ${sample}.${ref_fasta.baseName}.svision-pro.vcf.gz ${sample}.${ref_fasta.baseName}.svision-pro.bcf 
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
        .set { samples_ch }
    
    SVISIONPRO(samples_ch, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_REHEADER_SORT(SVISIONPRO.out, file(params.ref_fasta))
}
