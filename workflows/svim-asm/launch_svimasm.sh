#!/bin/bash
#SBATCH --job-name=launch_svimasm
#SBATCH --output=launch_svimasm-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=28-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
REF_FASTA="GRCh38_no_alt.fa"
REF_KMERS="repetitive_k19.txt"
OUTDIR="results"

nextflow run main.nf \
	-ansi-log false \
	-profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
	--ref_fasta ${REF_FASTA} \
	--ref_kmers ${REF_KMERS} \
	--outdir ${OUTDIR}
