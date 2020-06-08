#!/usr/bin/env python3

"""
Take bed-like files and present the data per base pair (and not by range)
Alexandre Pellan Cheng 2018
"""
from pyfaidx import Fasta
import os
import re
import sys

def one_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth):
    if hg19_seq == 'C' and hg19_seq_next == 'G': #most likely
        output.write(chr+"\t"+str(start)+"\t"+str(start+1)+"\t"+pmeth)

    elif hg19_seq_prev == 'C' and hg19_seq == 'G': #less likely, but possible (reverse strand or weird formatting)
        output.write(chr+"\t"+str(start-1)+"\t"+str(start)+"\t"+pmeth)

def two_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth):
    if hg19_seq == 'CG': #most likely
        output.write(chr+"\t"+str(start)+"\t"+str(start+1)+"\t"+pmeth)

    elif hg19_seq_prev == 'CG':
        output.write(chr+"\t"+str(start-1)+"\t"+str(start)+"\t"+pmeth)

def three_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth):
    prevbp = hg19_seq_prev[0]
    bp1 = hg19_seq[0]
    bp2 = hg19_seq[1]
    bp3 = hg19_seq[2]
    nextbp = hg19_seq_next[2]

    if bp1 == "C" and bp2 == "G":
        output.write(chr+"\t"+str(start)+"\t"+str(start+1)+"\t"+pmeth)
    elif bp2 == "C" and bp3 == "G":
        output.write(chr+"\t"+str(start+1)+"\t"+str(start+2)+"\t"+pmeth)
    elif bp3 == "C" and nextbp == "G":
        output.write(chr+"\t"+str(start+2)+"\t"+str(start+3)+"\t"+pmeth)
    elif prevbp == 'C' and bp1 == 'G':
        output.write(chr+"\t"+str(start-1)+"\t"+str(start)+"\t"+pmeth)

def four_or_more_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth):
    if bool(re.match('^(?!.*(CC|GG))[CG]*$', hg19_seq)):
        for bp in hg19_seq:
            if bp=="C":
                output.write(chr+"\t"+str(start)+"\t"+str(start+1)+"\t"+pmeth)
            start=start+1

def scaleback(tissuefile, chr_lengths_path, hg19_path, outfile):
    hg19 = Fasta(hg19_path)
    chr_len={}
    with open(chr_lengths_path) as f:
        for line in f:
            (key, val) = line.split()
            chr_len[key] = val
    in_file=tissuefile
    with open(in_file, "r") as referencefile:
        first_line_header=False
        if "bed" in referencefile.readline():
            first_line_header=True
    with open(in_file, "r") as referencefile, open(outfile, "w") as output:
        if first_line_header:
            next(referencefile)
        for line in referencefile:
            chr, start, end, pmeth = line.split("\t")
            start = int(float(start))
            end = int(float(end))
            if chr in chr_len:
                if end < int(chr_len[chr]) and start > 0:
                    hg19_seq = str(hg19[chr][start:end]).upper()
                    hg19_seq_prev = str(hg19[chr][(start-1):(end-1)]).upper()
                    hg19_seq_next = str(hg19[chr][(start+1):(end+1)]).upper()

                    if len(hg19_seq) == 1:
                        one_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth)

                    elif len(hg19_seq) == 2:
                        two_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth)

                    elif len(hg19_seq) == 3:
                        three_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth)

                    elif len(hg19_seq) > 3: #Only handles cases or repetitive CGCGCG... TEST THIS SECTION
                        four_or_more_nuc(hg19_seq, hg19_seq_prev, hg19_seq_next, output, chr, start, end, pmeth)

def main():
    tissuefile=sys.argv[1]
    chr_len=sys.argv[2]
    hg19_path=sys.argv[3]
    outfile = sys.argv[4]
    scaleback(tissuefile, chr_len, hg19_path, outfile)

if __name__ == '__main__':
    main()
