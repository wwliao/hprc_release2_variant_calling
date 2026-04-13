#!/bin/bash
#SBATCH --job-name=launch_kanpig
#SBATCH --output=launch_kanpig-%j.out
#SBATCH --partition=pi_hall
#SBATCH --constraint=nogpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=1-00:00:00

module purge
module load Nextflow/24.04.4

SAMPLE_SHEET="samplesheet.csv"
OUTDIR="results"
REF_FASTA="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa"
REF_FAI="/gpfs/gibbs/pi/ycgh/wl474/resources/reference_genomes/GRCh38_no_alt/GRCh38_no_alt.fa.fai"
PLOIDY_MALE_BED="/gpfs/gibbs/pi/ycgh/wl474/projects/hprc_r2/variant_benchmarking/230_samples/svdss/GRCh38_no_alt/kanpig/ploidy_beds/grch38_male.bed"
PLOIDY_FEMALE_BED="/gpfs/gibbs/pi/ycgh/wl474/projects/hprc_r2/variant_benchmarking/230_samples/svdss/GRCh38_no_alt/kanpig/ploidy_beds/grch38_female.bed"

nextflow run main.nf \
    -ansi-log true \
    -profile mccleary \
    --sample_sheet ${SAMPLE_SHEET} \
    --ref_fasta ${REF_FASTA} \
    --ref_fai ${REF_FAI} \
    --ploidy_male_bed ${PLOIDY_MALE_BED} \
    --ploidy_female_bed ${PLOIDY_FEMALE_BED} \
    --outdir ${OUTDIR}
