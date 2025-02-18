#!/usr/bin/env nextflow

process SAMTOOLS_MERGE_SORT {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam1), path(bam2)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.sorted.bam"), path("${sample}.${ref_fasta.baseName}.sorted.bam.bai")
    
    script:
    """
    samtools merge -@ ${task.cpus} - ${bam1} ${bam2} | samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sample}.${ref_fasta.baseName}.sorted.bam##idx##${sample}.${ref_fasta.baseName}.sorted.bam.bai -
    """
}

process SVISIONPROASM {
    tag "${sample}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    
    output:
    tuple val(sample), path("${sample}.svision_pro_v2.4.s1.vcf")

    script:
    """
    SVision-pro --process_num ${task.cpus} --preset asm --min_supp 1 --min_mapq 5 --min_sv_size 30 --detect_mode germline --img_size 512 --model_path /opt/svision-pro/src/pre_process/model_liteunet_256_8_16_32_32_32.pth --sample_name ${sample} --target_path ${bam} --genome_path ${ref_fasta} --out_path \${PWD}
    """
}

process BCFTOOLS_REHEADER_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.svision-pro-asm.vcf.gz"), path("${sample}.${ref_fasta.baseName}.svision-pro-asm.vcf.gz.tbi")

    script:
    """
    echo "${sample}" > sample.txt
    bcftools reheader -s sample.txt -o ${sample}.${ref_fasta.baseName}.svision-pro-asm.bcf ${vcf}
    bcftools sort --write-index=tbi -Oz -o ${sample}.${ref_fasta.baseName}.svision-pro-asm.vcf.gz ${sample}.${ref_fasta.baseName}.svision-pro-asm.bcf 
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.bam1), file(row.bam2)) }
        .set { samples_ch }

    SAMTOOLS_MERGE_SORT(samples_ch, file(params.ref_fasta))
    
    SVISIONPROASM(SAMTOOLS_MERGE_SORT.out, file(params.ref_fasta), file(params.ref_fai))

    BCFTOOLS_REHEADER_SORT(SVISIONPROASM.out, file(params.ref_fasta))
}
