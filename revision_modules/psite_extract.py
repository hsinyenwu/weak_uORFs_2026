#!/usr/bin/env python3
"""P-site extraction from SAM text on stdin (samtools view -q 255 [-L roi.bed] file.bam).
usage: psite_extract.py RULE > counts.tsv   (columns: chr, pos(1-based), strand, count)
Offsets are applied along the read from its 5' end, walking across spliced blocks (N) and deletions (D)."""
import sys, re
RULES = {
    "A": {24: 8, 25: 9, 26: 10, 27: 11, 28: 12},                                   # RiboTaper settings used by the lab (24-28 nt)
    "B": {23: 7, 24: 8, 25: 9, 26: 10, 27: 11, 28: 12, 29: 12, 30: 13},           # RiboseQC cutoffs for the merged libraries
    "C": {L: L - 16 for L in range(21, 37)},                                       # length-dependent, 3' anchored
    "D": {L: 12 for L in range(20, 37)},                                           # fixed 12
}
off = RULES[sys.argv[1]]
cig_re = re.compile(r"(\d+)([MIDNSHP=X])")
counts = {}
n_in = n_used = 0
for line in sys.stdin:
    n_in += 1
    f = line.split("\t", 10)
    L = len(f[9])
    o = off.get(L)
    if o is None:
        continue
    flag = int(f[1]); chrom = f[2]; ref = int(f[3]); cigar = f[5]
    blocks = []
    for n, op in cig_re.findall(cigar):
        n = int(n)
        if op in "M=X":
            blocks.append((ref, ref + n - 1)); ref += n
        elif op in "DN":
            ref += n
        # I, S, H, P: no reference consumption
    if not blocks:
        continue
    if flag & 16:
        strand = "-"; p = None
        for s, e in reversed(blocks):
            w = e - s + 1
            if o < w:
                p = e - o; break
            o -= w
    else:
        strand = "+"; p = None
        for s, e in blocks:
            w = e - s + 1
            if o < w:
                p = s + o; break
            o -= w
    if p is None:
        continue
    key = (chrom, strand, p)
    counts[key] = counts.get(key, 0) + 1
    n_used += 1
out = sys.stdout
for (chrom, strand, p), c in counts.items():
    out.write(f"{chrom}\t{p}\t{strand}\t{c}\n")
sys.stderr.write(f"reads in: {n_in}  used: {n_used}  positions: {len(counts)}\n")
