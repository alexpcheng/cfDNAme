#!/usr/bin/env python3

"""
Alexandre Pellan Cheng 2019
"""

import argparse
import glob
import os
from joblib import Parallel, delayed

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=str)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--metilene_output_path')
    parser.add_argument('--outdir')
    args = parser.parse_args()
    out = args.out
    metilene_path = args.metilene_output_path
    outdir = args.outdir
    threads = args.threads
    return(out, threads)


def tissue_DMR_per_chr(chr, tissue_group, metilene_path):
    files_to_compare = glob.glob(f'{metilene_path}/{chr}/*{tissue_group}_*')
    multi_intersect = 'multiIntersectBed -i ' + ' '.join(files_to_compare)
    awk_command = "awk -v OFS='\t' '{ for(i=6; i<=NF;i++) j+=$i; print $0,j; j=0 }'"
    os.system(multi_intersect + ' | ' + awk_command + f' > {tissue_group}.intersect.{chr}.tmp')

def merge_DMR_per_chr(autosomal_chr, tissue_group, outdir):
    out_file = f'{outdir}/{tissue_group}_multi_intersect.bed'
    files_to_cat = ' '.join([f'{tissue_group}.intersect.'+chr+'.tmp' for chr in sorted(autosomal_chr)])
    cat_command = f'cat {files_to_cat} > {out_file}'
    os.system(cat_command)
    os.system(f'rm {tissue_group}.intersect.*.tmp')

def DMR_in_most_groups(min_groups, tissue_group):
    new_out = f'sandbox/MethylMatrix/golden_markers/{tissue_group}_multi_intersect.filt.bed'
    with open(f'sandbox/MethylMatrix/golden_markers/{tissue_group}_multi_intersect.bed') as f, open(new_out, 'w') as w:
        for line in f:
            if int(line.strip().split('\t')[25]) >= min_groups:
                chr, start, end = line.split('\t')[0:3]
                num_groups = line.split('\t')[-1]
                w.write('\t'.join([chr, start, end+'\n']))
    os.system(f'bedtools merge -i {new_out} > {tissue_group}.merged.tmp && mv {tissue_group}.merged.tmp {new_out}')

def merge_all_golden_markers(out):
    command = f'cat sandbox/MethylMatrix/golden_markers/*multi_intersect.filt.bed | bedtools sort -i stdin | bedtools merge -i - > {out}'
    os.system(command)

def main():
    out, threads = parse_arguments()
    autosomal_chr = ['chr' + str(i) for i in range(1,23)]
    tissue_groups = ['G' + str(i) for i in range(1,22)]
    min_groups = 15
    print('tissue_DMR_per_chr')
    for tissue_group in tissue_groups:
        Parallel(n_jobs=threads)(delayed(tissue_DMR_per_chr)(chr, tissue_group) for chr in autosomal_chr)
    print('merge_DMR_per_chr')
    Parallel(n_jobs=threads)(delayed(merge_DMR_per_chr)(autosomal_chr, tissue_group) for tissue_group in tissue_groups)
    print('DMR_in_most_groups')
    Parallel(n_jobs=threads)(delayed(DMR_in_most_groups)(min_groups, tissue_group) for tissue_group in tissue_groups)
    print('merge_all_golden_markers')
    print(out)
    merge_all_golden_markers(out)
if __name__ == '__main__':
    main()
