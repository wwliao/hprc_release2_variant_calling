# Variant Calling for the HPRC Release 2

This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 231 individuals included in the Human Pangenome Reference Consortium (HPRC) Release 2.

## Overview

We used **Winnowmap (v2.03)** to align assemblies and HiFi reads to two reference genomes: **GRCh38\_no\_alt** and **CHM13v2**. However, for **PAV**, we used **Minimap2 (v2.26)** instead of Winnowmap. For HiFi reads with CpG methylation information (MM and ML tags), we retained the methylation data in the aligned BAM files.

Variants were then called using various tools. For structural variant (SV) callers, we relaxed the calling criteria to maximize recall by setting the following parameters (where applicable):

- Minimum MAPQ: 5
- Minimum read support: 3
- Minimum SV length: 30 bp

The table below summarizes the current status for each variant caller:

| Caller                 | Method Type    | GRCh38\_no\_alt | CHM13v2 |
| ---------------------: | :------------: | :-------------: | :-----: |
| CuteSV-asm (v2.1.1)    | Assembly-based | 231/231         | 231/231 |
| Dipcall (v0.3)         | Assembly-based | 231/231         | 231/231 |
| PAV (v2.4.6)           | Assembly-based | 231/231         | 231/231 |
| SVIM-asm (v1.0.3)      | Assembly-based | 231/231         | 231/231 |
| SVision-pro-asm (v2.4) | Assembly-based | 231/231         | 231/231 |
| DeepVariant (v1.6.1)   | HiFi-based     | 231/231         | 231/231 |
| CuteSV (v2.1.1)        | HiFi-based     | 231/231         | 231/231 |
| DeBreak (v1.3)         | HiFi-based     | 231/231         | 231/231 |
| Delly (v1.3.2)         | HiFi-based     | 231/231         | 231/231 |
| PBSV (v2.10.0)         | HiFi-based     | 231/231         | 231/231 |
| Sawfish (v0.12.8)      | HiFi-based     | 231/231         | 231/231 |
| Sniffles (v2.5.3)      | HiFi-based     | 231/231         | 231/231 |
| SVDSS (v2.0.0)         | HiFi-based     | 231/231         | 231/231 |
| SVIM (v2.0.0)          | HiFi-based     | 231/231         | 231/231 |
| SVision-pro (v2.4)     | HiFi-based     | 231/231         | 231/231 |

### Notes

- We modified `diploid_calling.py` in **CuteSV-asm (v2.1.1)** to support user-defined haplotype names. The modified version is available [here](https://github.com/wwliao/cuteSV).
- The GT field in **SVision-pro-asm (v2.4)** VCF output is incorrect ([see issue](https://github.com/songbowang125/SVision-pro/issues/15)). If needed, use the RNAMES field in INFO to infer the genotype.
- **DeepVariant (v1.8.0)** had a bug when used with CHM13v2 ([see issue](https://github.com/google/deepvariant/issues/912#issuecomment-2552635974)). To ensure consistency, we switched to **v1.6.1** for both reference genomes as a workaround.
- **CuteSV (v2.1.1)** occasionally reports variants with a position of zero ([see issue](https://github.com/tjiangHIT/cuteSV/issues/147)). We applied a post-processing step to remove these records before sorting to avoid issues with BCFtools.
- **SVDSS (v2.0.0)** currently calls only INS and DEL. The GT field in the VCF output is unreliable and should not be used.
- We modified **SVIM (v2.0.0)** to fix errors:

  1. Replaced `scipy`'s `linkage` with `fastcluster`'s `linkage` for hierarchical clustering.
  2. Updated `legendHandles` to `legend_handles` for Matplotlib compatibility.

  The modified version can be found [here](https://github.com/wwliao/svim).

## Index Files

| Index Type                       | Description                                            | File Name |
| -------------------------------: | :----------------------------------------------------- | :-------- |
| Assemblies                       | List of all assemblies included in HPRC Release 2      | [assemblies_pre_release_v0.6.1.index.csv](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/assemblies_pre_release_v0.6.1.index.csv) |
| HiFi Reads                       | List of all PacBio HiFi reads used in variant calling  | [hifi_reads.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/data_tables/hifi_reads.index.csv)|
| Assembly-to-Reference Alignments | List of assembly alignments to reference genomes       | [assembly_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/data_tables/assembly_alignments.index.csv)|
| HiFi-to-Reference Alignments     | List of HiFi read alignments to reference genomes      | [hifi_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/data_tables/hifi_alignments.index.csv)|
| Variant Callsets                 | List of all variant callsets generated for each sample | [variant_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/data_tables/variant_callsets.index.csv)|

### How to Download Files

To download files from the AWS S3 bucket, install the [AWS Command Line Interface (AWS CLI)](https://aws.amazon.com/cli/) if you haven't already, then run:

```bash
aws s3 --no-sign-request cp <s3_path> .
```

## Reference Genomes

We provide two reference genome packages: [`GRCh38_no_alt.tar.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/references/GRCh38_no_alt.tar.gz) and [`CHM13v2.tar.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/references/CHM13v2.tar.gz). Each package contains all necessary files for variant calling workflows.

### File Structure

After extracting `GRCh38_no_alt.tar.gz` or `CHM13v2.tar.gz`, you will find the following files:

- **`<reference>.fa`**: The FASTA file containing the reference genome.
- **`<reference>.fa.fai`**: Index file for the FASTA reference genome.
- **`<reference>.fmd`**: FMD index file for the FASTA reference genome.
- **`<reference>.PAR.bed`**: BED file specifying pseudo-autosomal regions (PARs).
- **`<reference>.expected_cn.XX.bed`**: BED file for expected copy numbers in female samples.
- **`<reference>.expected_cn.XY.bed`**: BED file for expected copy numbers in male samples.
- **`<reference>.TRF.bed`**: BED file annotating tandem repeats.
- **`repetitive_k15.txt`**: Text file listing repetitive k-mers (k=15) pre-computed using [meryl](https://github.com/marbl/meryl).
- **`repetitive_k19.txt`**: Text file listing repetitive k-mers (k=19) pre-computed using [meryl](https://github.com/marbl/meryl).

Replace `<reference>` with `GRCh38_no_alt` or `CHM13v2`, depending on the genome package you are using.

### GRCh38\_no\_alt.fa

The `GRCh38_no_alt.fa` file is derived from [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). It was decompressed and renamed to `GRCh38_no_alt.fa`. This version:

- Excludes ALT contigs.
- Has the PARs on chrY hard-masked.
- Uses the rCRS mitochondrial sequence.
- Includes the Epstein-Barr Virus (EBV) sequence.

For more details, see Heng Li’s blog post: [_Which human reference genome to use?_](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).

### CHM13v2.fa

The `CHM13v2.fa` file is derived from [`chm13v2.0_maskedY_rCRS.fa.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz). It was decompressed, modified to include the EBV sequence, and renamed to `CHM13v2.fa`. This version:

- Has the PARs on chrY hard-masked.
- Uses the rCRS mitochondrial sequence.
- Includes the EBV sequence.

