# Variant Calling for the HPRC Release 2

This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 231 individuals in the Human Pangenome Reference Consortium Release 2 (HPRC R2). It includes workflows for aligning assemblies and HiFi reads to reference genomes, as well as workflows for calling small and structural variants from these alignments using multiple tools. Additionally, [index files](#index-files) are provided for easy access to assemblies, HiFi reads, alignments, variant callsets, and merged callsets.

## Overview

We used **Winnowmap (v2.03)** to align assemblies and HiFi reads to two reference genomes: **GRCh38\_no\_alt** and **CHM13v2**. However, for **PAV**, we used **Minimap2 (v2.26)** instead of Winnowmap. For HiFi reads with CpG methylation information (MM and ML tags), we retained the methylation data in the aligned BAM files.

Variants were then called using various tools. For structural variant (SV) callers, we relaxed the calling criteria to maximize recall by setting the following parameters (where applicable):

- Minimum MAPQ: 5
- Minimum read support: 3
- Minimum SV length: 30 bp

In addition to the 231 HPRC R2 individuals, we also generated variant calls for an HPRC R2–equivalent version of HG002. The HG002 sequencing data were downsampled to match the typical coverage of HPRC R2 samples and processed through the same automated assembly and polishing pipeline. Variant calling was then performed on both the resulting assembly and the downsampled HiFi reads using the same workflows as for the other samples. This callset can be used for benchmarking against the GIAB T2T-HG002 Q100 v1.1 truth set.

## Variant Callers

Variant callers are grouped by input type. Most callers report SVs, some are joint callers that report both small variants and SVs, and DeepVariant reports only small variants. cuteSV can operate on both assemblies and HiFi reads; we use the suffix "-asm" for assembly-based runs.

### Assembly-based callers

| Caller          | Version | Type  | Notes                              |
|----------------:|:-------:|:-----:|:-----------------------------------|
| dipcall         | v0.3    | joint |                                    |
| PAV             | v2.4.6  | joint | Uses Minimap2 instead of Winnowmap |
| cuteSV-asm      | v2.1.1  | SV    | Modified implementation            |
| SVIM-asm        | v1.0.3  | SV    |                                    |

### HiFi-based callers

| Caller        | Version | Type  | Notes                             |
|--------------:|:-------:|:-----:|:----------------------------------|
| longcallD     | v0.0.4  | joint | Phased BAMs (HP/PS tags)          |
| DeepVariant   | v1.6.1  | small | v1.8.0 incompatible with CHM13v2  |
| cuteSV        | v2.1.1  | SV    | Post-processing applied           |
| DeBreak       | v1.3    | SV    |                                   |
| DELLY         | v1.3.2  | SV    |                                   |
| pbsv          | v2.10.0 | SV    |                                   |
| sawfish       | v0.12.8 | SV    |                                   |
| Sniffles2     | v2.5.3  | SV    |                                   |
| SVDSS         | v2.0.0  | SV    | Regenotyped                       |
| SVIM          | v2.0.0  | SV    | Modified implementation           |

### Notes

- **cuteSV-asm**: We modified `diploid_calling.py` to support user-defined haplotype names. The modified version is available [here](https://github.com/wwliao/cuteSV).

- **longcallD**: We used `--out-bam` to generate phased BAM files (HP and PS tags). Since the original input BAMs (listed in [hifi_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_alignments.index.csv)) are already available in the S3 bucket, we stored only extracted HP/PS tags in TSV format to reduce storage usage. Phased BAMs can be reconstructed using [`restore_phased_bam.py`](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/workflows/longcalld/restore_phased_bam.py) (requires `pysam`).

- **DeepVariant**: Version 1.8.0 has a bug when used with CHM13v2 ([see issue](https://github.com/google/deepvariant/issues/912#issuecomment-2552635974)), so we used v1.6.1 for both reference genomes.

- **cuteSV**: Occasionally reports variants with position zero ([see issue](https://github.com/tjiangHIT/cuteSV/issues/147)). These records were removed before sorting to avoid issues with BCFtools.

- **SVDSS**: Reports only insertions and deletions. The GT field is unreliable, so we used [kanpig](https://github.com/ACEnglish/kanpig) for regenotyping.

- **SVIM**: Modified by replacing `scipy.cluster.hierarchy.linkage` with `fastcluster.linkage` and updating `legendHandles` to `legend_handles` for compatibility with newer Matplotlib versions. The modified version is available [here](https://github.com/wwliao/svim).

## Merged Callsets

In addition to per-caller variant callsets, we generated **merged callsets** for each sample by integrating results from multiple callers. These merged callsets provide a more comprehensive set of variants by combining evidence across methods. The detailed methodology for constructing merged callsets, including their use in benchmarking, is described in the [graph variant benchmarking repository](https://github.com/wwliao/hprc_release2_graph_variant_benchmarking).

Merged callsets were generated for *autosomal chromosomes only*, as many callers do not account for sex when calling variants on sex chromosomes, which can result in inaccurate genotypes. Variants are further restricted to confident genomic regions defined by dipcall BED files for each sample.

An index of merged callsets is provided in the [Index Files](#index-files) section.

**Note**: A merged callset was not generated for HG00272, as this sample was not included in the HPRC R2 pangenome graph due to a likely large-scale misassembly on chromosome X.

## Index Files

| Index Type                       | Description                                            | File Name |
| -------------------------------: | :----------------------------------------------------- | :-------- |
| Assemblies                       | List of all assemblies used in variant calling         | [assemblies.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/assemblies.index.csv) |
| HiFi Reads                       | List of all PacBio HiFi reads used in variant calling  | [hifi_reads.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_reads.index.csv)|
| Assembly-to-Reference Alignments | List of assembly alignments to reference genomes       | [assembly_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/assembly_alignments.index.csv)|
| HiFi-to-Reference Alignments     | List of HiFi read alignments to reference genomes      | [hifi_alignments.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/hifi_alignments.index.csv)|
| Variant Callsets                 | List of all variant callsets generated for each sample | [variant_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/variant_callsets.index.csv)|
| Merged Callsets                  | List of merged variant callsets for each sample        | [merged_callsets.index.csv](https://github.com/wwliao/hprc_release2_variant_calling/blob/main/index_files/merged_callsets.index.csv)|

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

