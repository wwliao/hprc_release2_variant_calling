# Variant Calling for the HPRC Release 2

This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 216 individuals included in the Human Pangenome Reference Consortium (HPRC) Release 2.

## Overview

We used **Winnowmap (v2.03)** to align assemblies and HiFi reads to two reference genomes: **GRCh38\_no\_alt** and **CHM13v2**. However, for **PAV**, we used **Minimap2 (v2.26)** instead of Winnowmap. For HiFi reads with CpG methylation information (MM and ML tags), we retained the methylation data in the aligned BAM files.

Variants were then called using various tools. For structural variant (SV) callers, we relaxed the calling criteria to maximize recall by setting the following parameters (where applicable):

- Minimum MAPQ: 5
- Minimum read support: 3
- Minimum SV length: 30 bp

The table below summarizes the current status for each variant caller:

| Caller               | GRCh38\_no\_alt | CHM13v2         |
| ---:                 | ---:            | ---:            |
| Dipcall (v0.3)       | 216/216         | 216/216         |
| PAV (v2.4.6)         | 216/216         | 216/216         |
| SVIM-asm (v1.0.3)    | 216/216         | 216/216         |
| DeepVariant (v1.6.1) | 216/216         | 216/216         |
| CuteSV (v2.1.1)      | 216/216         | 216/216         |
| DeBreak (v1.3)       | 216/216         | 216/216         |
| Delly (v1.3.2)       | 216/216         | 216/216         |
| PBSV (v2.10.0)       | 216/216         | 216/216         |
| Sawfish (v0.12.8)    | 216/216         | 216/216         |
| Sniffles (v2.5.3)    | 216/216         | 216/216         |
| SVIM (v2.0.0)        | 216/216         | 216/216         |

### Notes

- **DeepVariant (v1.8.0)** had a bug when used with CHM13v2 ([see issue](https://github.com/google/deepvariant/issues/912#issuecomment-2552635974)). To ensure consistency, we switched to **v1.6.1** for both reference genomes as a workaround.
- We modified **SVIM (v2.0.0)** to fix errors:

  1. Replaced `scipy`'s `linkage` with `fastcluster`'s `linkage` for hierarchical clustering.
  2. Updated `legendHandles` to `legend_handles` for Matplotlib compatibility.

  The modified version can be found in [this GitHub repository](https://github.com/wwliao/svim).

## Reference Genomes

We provide two reference genome packages: [`GRCh38_no_alt.zip`](https://drive.google.com/uc?id=10bh1CEv0ifHVv9nNHTXI0TlChiN7Bj6D) and [`CHM13v2.zip`](https://drive.google.com/uc?id=1XECM8XeWVLY3NZsYvCBM6ipU7GJ-trbW). Each package contains all necessary files for variant calling workflows.

### File Structure

After extracting `GRCh38_no_alt.zip` or `CHM13v2.zip`, you will find the following files:

- **`<reference>.fa`**: The FASTA file containing the reference genome.
- **`<reference>.fa.fai`**: Index file for the FASTA reference genome.
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

