#!/usr/bin/env python3
import argparse
import gzip

import pysam

parser = argparse.ArgumentParser(description="Restore HP/PS tags to a BAM file from a TSV file")
parser.add_argument("tsv_file", help="Input TSV file (can be .gz)")
parser.add_argument("input_bam", help="Original BAM file")
parser.add_argument("output_bam", help="Output BAM file with restored HP/PS tags")
args = parser.parse_args()

haplotype_info = {}
open_func = gzip.open if args.tsv_file.endswith(".gz") else open

with open_func(args.tsv_file, "rt") as f:
    next(f)
    for line in f:
        read_name, chrom, pos, hp, ps = line.strip().split("\t")
        pos = int(pos)
        hp = int(hp) if hp != "." else None
        ps = int(ps) if ps != "." else None

        haplotype_info[(read_name, chrom, pos)] = (hp, ps)

with pysam.AlignmentFile(args.input_bam, "rb") as in_bam, \
     pysam.AlignmentFile(args.output_bam, "wb", template=in_bam) as out_bam:

    for read in in_bam:
        # Ensure this is a primary alignment (not secondary and not supplementary)
        if not read.is_secondary and not read.is_supplementary:
            key = (read.query_name, read.reference_name, read.reference_start + 1)  # Convert 0-based BAM to 1-based
            if key in haplotype_info:
                hp, ps = haplotype_info[key]
                if hp is not None:
                    read.set_tag("HP", hp, value_type="i")
                if ps is not None:
                    read.set_tag("PS", ps, value_type="i")
                out_bam.write(read)
