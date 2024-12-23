#!/bin/bash
#SBATCH --job-name=launch_deepvariant
#SBATCH --output=launch_deepvariant-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=28-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
REF_FASTA="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa"
REF_FAI="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa.fai"
REF_PAR="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.PAR.bed"
OUTDIR="results"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --ref_fasta ${REF_FASTA} \
    --ref_fai ${REF_FAI} \
    --ref_par ${REF_PAR} \
    --outdir ${OUTDIR}
