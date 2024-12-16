#!/usr/bin/env nextflow

process WINNOWMAP_ASM {
    tag "${sample}_${hap_id}"

    input:
    tuple val(sample), val(hap_id), path(hap_fasta)
    path ref_fasta
    path ref_kmers

    output:
    tuple val(sample), val(hap_id), path("${sample}_${hap_id}.${ref_fasta.simpleName}.sam")

    script:
    """
    winnowmap -W ${ref_kmers} -a -x asm5 --cs -r 2k -t ${task.cpus} ${ref_fasta} ${hap_fasta} > ${sample}_${hap_id}.${ref_fasta.simpleName}.sam
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "${sample}_${hap_id}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), val(hap_id), path(sam)
    
    output:
    tuple val(sample), val(hap_id), path("${sam.baseName}.sorted.bam"), path("${sam.baseName}.sorted.bam.bai")
    
    script:
    """
    samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sam.baseName}.sorted.bam##idx##${sam.baseName}.sorted.bam.bai ${sam}
    """
}

process SVIMASM_DIPLOID {
    tag "${sample}"

    input:
    tuple val(sample), path(bam1), path(bai1), path(bam2), path(bai2)
    path ref_fasta

    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.svim-asm.vcf")

    script:
    """
    svim-asm diploid --interspersed_duplications_as_insertions --sample ${sample} ${sample}_svimasm_results ${bam1} ${bam2} ${ref_fasta}
    mv ${sample}_svimasm_results/variants.vcf ${sample}.${ref_fasta.simpleName}.svim-asm.vcf
    """
}

process BCFTOOLS_SORT_INDEX {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(vcf)

    output:
    tuple val(sample), path("${vcf}.gz"), path("${vcf}.gz.tbi")

    script:
    """
    bcftools sort --write-index=tbi -Oz -o ${vcf}.gz ${vcf}
    """
}

workflow {
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header:true)
        .map { row -> 
            [
                [sample: row.sample, hap_id: 'hap1', hap_fasta: file(row.hap1)],
                [sample: row.sample, hap_id: 'hap2', hap_fasta: file(row.hap2)]
            ]
        }
        .flatten()
        .map { it -> tuple(it.sample, it.hap_id, it.hap_fasta) }
        .set { samples_ch }

    // Align a haplotype to the reference genome using winnowmap
    WINNOWMAP_ASM(samples_ch, file(params.ref_fasta), file(params.ref_kmers))

    // Sort a SAM file into a BAM file and generate a BAI index file
    SAMTOOLS_SORT_INDEX(WINNOWMAP_ASM.out)

    SAMTOOLS_SORT_INDEX.out
        .map { sample, hap_id, bam, bai -> tuple(sample, tuple(hap_id, bam, bai)) }
        .groupTuple(by: 0, size: 2)
        .map { sample, hap_bams ->
            def sorted_hap_bams = hap_bams.sort { it[0] } // Sort by hap ID
            tuple(sample, sorted_hap_bams[0][1], sorted_hap_bams[0][2], sorted_hap_bams[1][1], sorted_hap_bams[1][2])
        }
        .set { grouped_bams_ch }

    // Call structural variants from diploid genome-genome alignments
    SVIMASM_DIPLOID(grouped_bams_ch, file(params.ref_fasta))

    // Sort a VCF file into a compressed VCF file and generate a TBI index file
    BCFTOOLS_SORT_INDEX(SVIMASM_DIPLOID.out)
}
