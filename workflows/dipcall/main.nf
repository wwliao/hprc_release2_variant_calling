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
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.dip.vcf.gz"), path("${sample}.${ref_fasta.baseName}.dipcall.bed"), emit: res
    
    script:
    if( sex == 'female' )
        """
        run-dipcall -t ${task.cpus} -W ${ref_kmers} ${sample}.${ref_fasta.baseName} ${ref_fasta} ${hap1} ${hap2} > ${sample}.${ref_fasta.baseName}.mak
        make -j 1 -f ${sample}.${ref_fasta.baseName}.mak
        mv ${sample}.${ref_fasta.baseName}.dip.bed ${sample}.${ref_fasta.baseName}.dipcall.bed
        """
    else if( sex == 'male' )
        """
        run-dipcall -t ${task.cpus} -x ${ref_par} -W ${ref_kmers} ${sample}.${ref_fasta.baseName} ${ref_fasta} ${hap1} ${hap2} > ${sample}.${ref_fasta.baseName}.mak
        make -j 1 -f ${sample}.${ref_fasta.baseName}.mak
        mv ${sample}.${ref_fasta.baseName}.dip.bed ${sample}.${ref_fasta.baseName}.dipcall.bed
        """
}

process SAMTOOLS_SORT {
    tag "${sample}"

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sam.baseName.replace('sam', 'sorted.bam')}"), path("${sam.baseName.replace('sam', 'sorted.bam.bai')}")

    script:
    """
    gunzip ${sam} | samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sam.baseName.replace('sam', 'sorted.bam')}##idx##${sam.baseName.replace('sam', 'sorted.bam.bai')} -
    """
}

process BCFTOOLS_REHEADER {
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(bed)

    output:
    tuple val(sample), path("${vcf.replace('dip', 'dipcall')}"), path(bed)

    script:
    """
    echo "syndip ${sample}" > sample_map.txt
    bcftools reheader -s sample_map.txt -o ${vcf.replace('dip', 'dipcall')} ${vcf}
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
    
    BCFTOOLS_REHEADER(DIPCALL.out.res)

    DIPCALL.out.sam
    .flatMap { sample, hap1_sam, hap2_sam -> 
        // Split hap1 and hap2 SAM files into separate tuples
        [tuple(sample, hap1_sam), tuple(sample, hap2_sam)]
    }
    .set { samfiles_ch }

    SAMTOOLS_SORT(samfiles_ch)
}
