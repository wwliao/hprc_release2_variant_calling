#!/usr/bin/env nextflow

process DIPCALL {
    tag "${sample}"
    
    input:
    tuple val(sample), val(sex), path(hap1), path(hap2)
    path ref_fasta
    path ref_fai
    path ref_kmers
    path ref_par
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.hap1.sam.gz"), path("${sample}.${ref_fasta.baseName}.hap2.sam.gz"), emit: sam
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.dip.vcf.gz"), path("${sample}.${ref_fasta.baseName}.dip.bed"), emit: dip
    
    script:
    def half_cpus = (task.cpus + 1).intdiv(2)
    if( sex == 'female' )
        """
        run-dipcall -t ${half_cpus} -W ${ref_kmers} ${sample}.${ref_fasta.baseName} ${ref_fasta} ${hap1} ${hap2} > ${sample}.${ref_fasta.baseName}.mak
        make -j2 -f ${sample}.${ref_fasta.baseName}.mak
        """
    else if( sex == 'male' )
        """
        run-dipcall -t ${half_cpus} -x ${ref_par} -W ${ref_kmers} ${sample}.${ref_fasta.baseName} ${ref_fasta} ${hap1} ${hap2} > ${sample}.${ref_fasta.baseName}.mak
        make -j2 -f ${sample}.${ref_fasta.baseName}.mak
        """
}

process BCFTOOLS_REHEADER_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(bed)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.dipcall.vcf.gz"), path("${sample}.${ref_fasta.baseName}.dipcall.vcf.gz.tbi"), path("${sample}.${ref_fasta.baseName}.dipcall.bed")

    script:
    """
    echo "${sample}" > sample.txt
    bcftools reheader -s sample.txt -o ${sample}.${ref_fasta.baseName}.dipcall.bcf ${vcf}
    bcftools sort --write-index=tbi -Oz -o ${sample}.${ref_fasta.baseName}.dipcall.vcf.gz ${sample}.${ref_fasta.baseName}.dipcall.bcf 
    mv ${bed} ${sample}.${ref_fasta.baseName}.dipcall.bed
    """
}

process SAMTOOLS_SORT {
    tag "${sample}_${hap}"

    input:
    tuple val(sample), val(hap), path(sam)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.${hap}.sorted.bam"), path("${sample}.${ref_fasta.baseName}.${hap}.sorted.bam.bai")

    script:
    """
    zcat ${sam} | samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sample}.${ref_fasta.baseName}.${hap}.sorted.bam##idx##${sample}.${ref_fasta.baseName}.${hap}.sorted.bam.bai -
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.sex, file(row.hap1), file(row.hap2)) }
        .set { samples_ch }

    DIPCALL(samples_ch, file(params.ref_fasta), file(params.ref_fai), file(params.ref_kmers), file(params.ref_par))
    
    BCFTOOLS_REHEADER_SORT(DIPCALL.out.dip, file(params.ref_fasta))

    DIPCALL.out.sam
    .flatMap { sample, hap1_sam, hap2_sam -> 
        // Split hap1 and hap2 SAM files into separate tuples
        [tuple(sample, 'hap1', hap1_sam), tuple(sample, 'hap2', hap2_sam)]
    }
    .set { samfiles_ch }

    SAMTOOLS_SORT(samfiles_ch, file(params.ref_fasta))
}
