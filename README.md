# Variant Calling for the HPRC Release 2
This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 216 individuals included in the Human Pangenome Reference Consortium (HPRC) Release 2.

## Overview
We used **Winnowmap (v2.03)** to align assemblies and HiFi reads to two reference genomes: **GRCh38\_no\_alt** and **CHM13v2**. However, for **PAV**, we used **minimap2 (v2.26)** instead of Winnowmap. Variants were then called using various tools. The table below summarizes the current status for each variant caller:


| Caller               | GRCh38\_no\_alt | CHM13v2         |
| ---:                 | ---:            | ---:            |
| Dipcall (v0.3)       | 216/216         | 216/216         |
| PAV (v2.4.6)         |   ?/216         |   ?/216         |
| SVIM-asm (v1.0.3)    | 216/216         | 216/216         |
| DeepVariant (v1.8.0) | 216/216         | Bug encountered |
| CuteSV (v2.1.1)      | 216/216         | 196/216         |
| DeBreak (v1.3)       |    -            |    -            |
| Delly (v1.3.2)       | 216/216         | 196/216         |
| PBSV (v2.10.0)       | 216/216         | 196/216         |
| Sawfish (v0.12.8)    | 216/216         | 196/216         |
| Sniffles (v2.5.3)    | 216/216         | 196/216         |
| SVIM (v2.0.0)        | 216/216         | 196/216         |

### Notes
- A bug in **PAV (v2.4.6)** was resolved ([see issue](https://github.com/EichlerLab/pav/issues/63#issuecomment-2510950978)), and the workflow has now been initiated.
- **DeepVariant (v1.8.0)** encountered a bug for CHM13v2 ([see issue](https://github.com/google/deepvariant/issues/912#issuecomment-2552635974)). We plan to switch to **v1.6.1** as a workaround.
- The **DeBreak (v1.3)** workflow has not been initiated for either reference genome due to high resource requirements and excessive generation of intermediate files.
- We modified **SVIM (v2.0.0)** to fix errors:
  1. Replaced `scipy`'s `linkage` with `fastcluster`'s `linkage` for hierarchical clustering.
  2. Updated `legendHandles` to `legend_handles` for Matplotlib compatibility.

  The modified version can be found in [my GitHub repository](https://github.com/wwliao/svim).
