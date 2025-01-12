# Variant Calling for the HPRC Release 2
This repository contains Nextflow workflows for assembly-based and HiFi-based variant calling for 216 individuals included in the Human Pangenome Reference Consortium (HPRC) Release 2.

## Overview
We used Winnowmap (v2.03) to align assemblies and HiFi reads to two reference genomes: GRCh38 and CHM13. Variants were then called using various tools, as detailed below. The following table summarizes the current status for each variant caller:

| Caller               | GRCh38  | CHM13             |
| ---:                 | ---:    | ---:              |
| Dipcall (v0.3)       | 216/216 | 216/216           |
| PAV (v2.4.6)         |   ?/216 |   ?/216           |
| SVIM-asm (v1.0.3)    | 216/216 | 216/216           |
| DeepVariant (v1.8.0) | 216/216 | bug               |
| CuteSV (v2.1.1)      | 216/216 | 160/216           |
| DeBreak (v1.3)       |    -    |    -              |
| Delly (v1.3.2)       | 216/216 |   ?/216           |
| PBSV (v2.10.0)       | 216/216 | 160/216           |
| Sawfish (v0.12.8)    | 216/216 |   ?/216           |
| Sniffles (v2.5.3)    | 216/216 | 160/216           |
| SVIM (v2.0.0)        | 216/216 | 160/216           |

