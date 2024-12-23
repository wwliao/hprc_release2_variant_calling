#!/usr/bin/env nextflow

process DEEPVARIANT {
    tag "${sample}"
    
    input:
    tuple val(sample), val(sex), path(bam), path(bai)
    path ref_fasta
    path ref_fai
    path ref_par
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.baseName}.deepvariant.vcf.gz"), path("${sample}.${ref_fasta.baseName}.deepvariant.vcf.gz.tbi"), path("${sample}.${ref_fasta.baseName}.deepvariant.g.vcf.gz"), path("${sample}.${ref_fasta.baseName}.deepvariant.g.vcf.gz.tbi")
    
    script:
    if( sex == 'female' )
        """
        mkdir -p "\${PWD}/intermediate_results"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --ref "${ref_fasta}" \
            --sample_name ${sample} \
            --reads "${bam}" \
            --output_vcf "${sample}.${ref_fasta.baseName}.deepvariant.vcf.gz" \
            --output_gvcf "${sample}.${ref_fasta.baseName}.deepvariant.g.vcf.gz" \
            --num_shards "${task.cpus}" \
            --intermediate_results_dir="\${PWD}/intermediate_results"
        """

    else if( sex == 'male' )
        """
        mkdir -p "\${PWD}/intermediate_results"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --ref "${ref_fasta}" \
            --sample_name ${sample} \
            --haploid_contigs="chrX,chrY" \
            --par_regions_bed="${ref_par}" \
            --reads "${bam}" \
            --output_vcf "${sample}.${ref_fasta.baseName}.deepvariant.vcf.gz" \
            --output_gvcf "${sample}.${ref_fasta.baseName}.deepvariant.g.vcf.gz" \
            --num_shards "${task.cpus}" \
            --intermediate_results_dir="\${PWD}/intermediate_results"
        """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.sex, file(row.bam), file(row.bai)) }
        .set { samples_ch }

    DEEPVARIANT(samples_ch, file(params.ref_fasta), file(params.ref_fai), file(params.ref_par))
}
