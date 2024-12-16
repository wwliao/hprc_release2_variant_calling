#!/usr/bin/env nextflow

process SAMTOOLS_FASTQ {
    tag "${bam}"
    
    input:
    tuple val(sample), path(bam)
    
    output:
    tuple val(sample), path("${bam.baseName}.fastq.gz")
    
    script:
    """
    samtools fastq -@ ${task.cpus} -0 ${bam.baseName}.fastq.gz ${bam}
    """
}

process WINNOWMAP_HIFI {
    tag "${fastq}"
    
    input:
    tuple val(sample), path(fastq)
    path ref_fasta
    path ref_kmers
    
    output:
    tuple val(sample), val(run_id), path("${sample}.${run_id}.${ref_fasta.simpleName}.bam")
    
    script:
    parts = fastq.name.split('\\.')
    if (parts.size() < 2) {
        error "Unable to extract run_id from filename: ${fastq.name}"
    }
    run_id = parts[1]
    rg_id = (run_id + "//CCS").md5().substring(0, 8)

    """
    winnowmap -W ${ref_kmers} -t ${task.cpus} -R "@RG\\tID:${rg_id}\\tPL:PACBIO\\tDS:READTYPE=CCS\\tPU:${run_id}\\tSM:${sample}" -x map-pb -a -Y -L --eqx --cs ${ref_fasta} ${fastq} | samtools sort -m 4G -@ ${task.cpus} -o ${sample}.${run_id}.${ref_fasta.simpleName}.bam
    """
}

// Module to sort and merge BAM files
process SAMTOOLS_MERGE_SORT {
    tag "${sample}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample), path(bams)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.sorted.bam"), path("${sample}.${ref_fasta.simpleName}.sorted.bam.bai")
    
    script:
    """
    samtools merge -@ ${task.cpus} - ${bams} | samtools sort -m 4G -@ ${task.cpus} --write-index -o ${sample}.${ref_fasta.simpleName}.sorted.bam##idx##${sample}.${ref_fasta.simpleName}.sorted.bam.bai -
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample, file(row.hifi)) }
        .set { samples_ch }
   
    // Separate BAM and FASTQ files
    samples_ch
        .branch {
            bam: it[1].name.endsWith('.bam')
            fastq: it[1].name.endsWith('.fastq.gz')
        }
        .set { branches_ch }

    // Convert BAM to FASTQ
    SAMTOOLS_FASTQ(branches_ch.bam)
    
    // Combine converted FASTQ and original FASTQ
    SAMTOOLS_FASTQ.out
        .mix(branches_ch.fastq)
        .toSortedList { a, b -> a[0] <=> b[0] }  // Sort by sample name
        .flatMap()
        .set { fastqs_ch }
   
    // Align all FASTQ files to reference using winnowmap
    WINNOWMAP_HIFI(fastqs_ch, file(params.ref_fasta), file(params.ref_kmers))
    
    // Group aligned BAM files by sample and sort by run ID
    WINNOWMAP_HIFI.out
        .map { sample, run_id, bam -> tuple(sample, tuple(run_id, bam)) }
        .groupTuple(by: 0)
        .map { sample, files -> 
            tuple(
                sample,
                files.sort { a, b -> a[0] <=> b[0] }.collect { it[1] }
            )
        }
        .set { grouped_bams_ch }
    
    // Merge and sort BAM files for each sample
    SAMTOOLS_MERGE_SORT(grouped_bams_ch, file(params.ref_fasta))
}
