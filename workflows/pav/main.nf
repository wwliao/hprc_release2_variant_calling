#!/usr/bin/env nextflow

process PAV {
    tag "${sample}"
    
    input:
    tuple val(sample), path(hap1), path(hap2)
    path ref_fasta
    
    output:
    tuple val(sample), path("${sample}.${ref_fasta.simpleName}.pav.vcf.gz"), path("${sample}.${ref_fasta.simpleName}.pav.vcf.gz.tbi")
    
    script:
    """
    cat << EOF > config.json
    {
      "reference": "${ref_fasta}"
    }
    EOF

    cat << EOF > assemblies.tsv
    NAME\tHAP1\tHAP2
    ${sample}\t${hap1}\t${hap2}
    EOF

    /opt/pav/files/docker/run --notemp -c ${task.cpus}
    mv ${sample}.vcf.gz ${sample}.${ref_fasta.simpleName}.pav.vcf.gz
    mv ${sample}.vcf.gz.tbi ${sample}.${ref_fasta.simpleName}.pav.vcf.gz.tbi
    """
}

workflow {
    // Read the sample sheet
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.hap1), file(row.hap2)) }
        .set { samples_ch }

    PAV(samples_ch, file(params.ref_fasta))
}
