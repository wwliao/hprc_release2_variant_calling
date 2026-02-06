# Variant Calling for the HPRC Release 2

This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 231 individuals in the Human Pangenome Reference Consortium (HPRC) Release 2. It includes workflows for aligning assemblies and HiFi reads to reference genomes, as well as workflows for calling small and structural variants from these alignments using multiple tools. Additionally, [index files](#index-files) are provided for easy access to assemblies, HiFi reads, alignments, and variant callsets.

## Overview

We used **Winnowmap (v2.03)** to align assemblies and HiFi reads to two reference genomes: **GRCh38\_no\_alt** and **CHM13v2**. However, for **PAV**, we used **Minimap2 (v2.26)** instead of Winnowmap. For HiFi reads with CpG methylation information (MM and ML tags), we retained the methylation data in the aligned BAM files.

Variants were then called using various tools. For structural variant (SV) callers, we relaxed the calling criteria to maximize recall by setting the following parameters (where applicable):

- Minimum MAPQ: 5
- Minimum read support: 3
- Minimum SV length: 30 bp

## Variant Callers

Variant callers are grouped by input type. Most callers report SVs, some are joint callers that report both small variants and SVs, and DeepVariant reports only small variants. cuteSV and SVision-pro can operate on both assemblies and HiFi reads, so we use the suffix "-asm" for assembly-based runs.

### Assembly-based callers

- **dipcall (v0.3, joint)**
- **PAV (v2.4.6, joint)**
- **cuteSV-asm (v2.1.1, SV)**

    We modified `diploid_calling.py` to support user-defined haplotype names. The modified version is available [here](https://github.com/wwliao/cuteSV).

- **SVIM-asm (v1.0.3, SV)**
- **SVision-pro-asm (v2.4, SV; excluded from downstream analysis)**

    The GT field in the VCF output is incorrect ([see issue](https://github.com/songbowang125/SVision-pro/issues/15)). The genotype can be inferred from the RNAMES field in INFO if needed. Because the VCF is not sequence-resolved, this caller was *excluded* from downstream analysis.

### HiFi-based callers

- **longcallD (v0.0.4, joint)**

    We used `--out-bam` to generate phased BAM files, which include HP and PS tags. Since the original input BAMs (listed in [hifi_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_alignments.index.csv)) are already available in the S3 bucket, we saved only the extracted HP and PS tag information as a TSV file to reduce storage usage. If you need the phased BAMs, you can easily reconstruct them using the input BAMs and the TSV file with the provided [`restore_phased_bam.py`](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/workflows/longcalld/restore_phased_bam.py) script (requires `pysam`).

- **DeepVariant (v1.6.1, small)**

    DeepVariant was used for small variant calling. Version 1.8.0 had a bug when used with CHM13v2 ([see issue](https://github.com/google/deepvariant/issues/912#issuecomment-2552635974)), so we used version 1.6.1 for both reference genomes to maintain consistency.

- **cuteSV (v2.1.1, SV)**

    cuteSV occasionally reports variants with a position of zero ([see issue](https://github.com/tjiangHIT/cuteSV/issues/147)). We applied a post-processing step to remove these records before sorting to avoid issues with BCFtools.

- **DeBreak (v1.3, SV)**
- **DELLY (v1.3.2, SV)**
- **pbsv (v2.10.0, SV)**
- **sawfish (v0.12.8, SV)**
- **Sniffles2 (v2.5.3, SV)**
- **SVDSS (v2.0.0, SV)**

    SVDSS currently reports only insertions and deletions. The GT field in the VCF output is unreliable, so we used [`kanpig`](https://github.com/ACEnglish/kanpig) to regenotype these variants.

- **SVIM (v2.0.0, SV)**

    We modified SVIM to fix errors by replacing `scipy.cluster.hierarchy.linkage` with `fastcluster.linkage` and updating `legendHandles` to `legend_handles` for compatibility with newer Matplotlib versions. The modified version is available [here](https://github.com/wwliao/svim).

- **SVision-pro (v2.4, SV; excluded from downstream analysis)**

    The VCF output is not sequence-resolved, so this caller was *excluded* from downstream analysis.

## Index Files

| Index Type                       | Description                                            | File Name |
| -------------------------------: | :----------------------------------------------------- | :-------- |
| Assemblies                       | List of all assemblies included in HPRC Release 2      | [assemblies.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/assemblies.index.csv) |
| HiFi Reads                       | List of all PacBio HiFi reads used in variant calling  | [hifi_reads.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_reads.index.csv)|
| Assembly-to-Reference Alignments | List of assembly alignments to reference genomes       | [assembly_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/assembly_alignments.index.csv)|
| HiFi-to-Reference Alignments     | List of HiFi read alignments to reference genomes      | [hifi_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_alignments.index.csv)|
| Variant Callsets                 | List of all variant callsets generated for each sample | [variant_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/variant_callsets.index.csv)|

### How to Download Files

To download files from the AWS S3 bucket, install the [AWS Command Line Interface (AWS CLI)](https://aws.amazon.com/cli/) if you haven't already, then run:

```bash
aws s3 --no-sign-request cp <s3_path> .
```

## Reference Genomes

We provide two reference genome packages: [`GRCh38_no_alt.tar.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/references/GRCh38_no_alt.tar.gz) and [`CHM13v2.tar.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/40399FDD-59DE-43D1-B3A3-DFF0C6E64FAC--YALE_VARIANT_CALLS_R2/references/CHM13v2.tar.gz). Each package contains all necessary files for the workflows in this repository.

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

