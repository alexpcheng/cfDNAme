#!/usr/bin/env python3

"""
Filters rare SNPs and underrepresented blocks
Eilon Sharon 2017
"""

import os
import argparse
import gzip

import pandas as pd

def filter_snps(geno_df, low_freq_threshold, is_threshold_for_any_pop, min_snps_per_block, min_reads_per_block):
    """
    filterring rare SNPs and under-representted blocks 
    """
    
    geno_df.rename(index=str, columns={"AF": "G1K_AF"}, inplace=True)
    af_cols = [colname for colname in geno_df.columns if "AF" in colname]
    
    if is_threshold_for_any_pop:
        geno_df['not_rare_snps_tf'] = (geno_df.ix[:,af_cols].min(axis=1) >=low_freq_threshold)
    else:
        geno_df['not_rare_snps_tf'] =  geno_df.ix[:,'G1K_AF']>=low_freq_threshold
        
    # TODO hack for simulations (if not selected for simulation define as rare SNPs that are excluded from the simulations)
    if 'selected_sim_SNP' in geno_df.columns:
        geno_df['not_rare_snps_tf'] = (geno_df['not_rare_snps_tf'] == True) & (geno_df['selected_sim_SNP'] == True)
        
    geno_df_blocks = geno_df.groupby('genetic_blocks')

    geno_df['block_filter']  = False

    filterred_blocks = 0
    cum_dist = 0
    for block_i, block_data in geno_df_blocks:
        cur_num_snps =  block_data.not_rare_snps_tf.sum()
        cur_num_reads = block_data.g1000_Allele1_seq_cnt[block_data.not_rare_snps_tf].sum() + block_data.g1000_Allele2_seq_cnt[block_data.not_rare_snps_tf].sum()
        cur_pass_block_filter = cur_num_snps >= min_snps_per_block and cur_num_reads >= min_reads_per_block
        cur_dist_from_prev_block = block_data.genetic_avg_distance_FromPrevBlock.mean()
        
        if cur_pass_block_filter:
            geno_df.ix[geno_df['genetic_blocks'] == block_i,'block_filter'] = True
            geno_df.ix[geno_df['genetic_blocks'] == block_i,'genetic_avg_distance_FromPrevBlock'] = cum_dist + cur_dist_from_prev_block
            cum_dist = 0
        else:
            cum_dist = cum_dist +  cur_dist_from_prev_block
            filterred_blocks = filterred_blocks +1
    
    print("Filterred %d blocks and %d SNPs" % (filterred_blocks, (~ (geno_df['not_rare_snps_tf'] & geno_df['block_filter'])).sum() ))
    
    geno_df = geno_df.ix[geno_df['not_rare_snps_tf'] & geno_df['block_filter'] ,]
    
    return(geno_df)
    
    
    


###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Filtering SNPs and SNPs block (for example: rare SNPs, under-represented blocks)")

    # input files

    parser.add_argument("input_data_filename", type=str, help= "tab-separated values (tsv) file produced by the cfDNA1G pipeline. Each SNP is a raw.")

    #output files
    
    parser.add_argument("output_filename", type=str, 
                      help="tab-separated values (tsv) file. A filttered version of the input")
    
    
    
    low_freq_threshold = 1e-5
    is_threshold_for_any_pop = False
    min_snps_per_block = 100
    min_reads_per_block = 100
    
    
    # options

    parser.add_argument("-f", "--low_freq_threshold", dest='low_freq_threshold', default=1e-5, type=float,
                      help='minimal SNP frequency')
    
    parser.add_argument("-a", "--is_threshold_for_any_pop", dest='is_threshold_for_any_pop',  action='store_true', default=False,
                        help='The frequency filter will be applied to any population (instead of the total 1000 genomes). Default is False (no flag)')
    
    parser.add_argument("-g", "--min_gc_score", dest='min_gc_score', default=0.7, type=float,
                      help='minimal GCScore to consider (default 0.7)')
        
    parser.add_argument("-s", "--min_snps_per_block", dest='min_snps_per_block', default=100, type=int,
                      help='minimal number of SNPs per block (default 100)')

    parser.add_argument("-r", "--min_reads_per_block", dest='min_reads_per_block', default=100, type=int,
                      help='minimal number of reads per block (default 100')
    
    args = parser.parse_args()
    
    
    
    assert(args.low_freq_threshold >= 0.0 and args.low_freq_threshold <= 1.0)


    # running the function
    
    print('Loading file: %s' % (args.input_data_filename))
    geno_df = pd.read_table(args.input_data_filename, '\t')
    
    if args.min_gc_score > 0 and 'GCScore' in geno_df.columns:
        print('filtering SNPs with illumina GC score lower than %f' % (args.min_gc_score))
        geno_df = geno_df.ix[ geno_df['GCScore'] >= args.min_gc_score,]

    #print('Filtering any duplicated lines and removing duplicated SNPs (and g1000_ID,GCScore columns)')
    #TODO - hack because of duplictaed SNPs in final report

    #if 'g1000_ID' in geno_df.columns:
    #    geno_df.drop('g1000_ID', axis=1, inplace=True)

    #if 'GCScore' in geno_df.columns:
    #    geno_df.drop('GCScore', axis=1, inplace=True)

    #if 'genetic_dist_log' in geno_df.columns:
    #    geno_df.drop('genetic_dist_log', axis=1, inplace=True)

    print('Filtering any duplicated lines - removing duplicated SNPs ')
    #TODO - hack because of duplictaed SNPs in final report
    snp_identifying_cols = ['Chr', 'Position', 'Strand', 'g1000_Allele1', 'g1000_Allele2']
    geno_df.drop_duplicates(inplace=True, subset = snp_identifying_cols, keep='first')

    print('filtering SNPs and blocks')
    geno_df = filter_snps(geno_df, args.low_freq_threshold, args.is_threshold_for_any_pop, args.min_snps_per_block, args.min_reads_per_block)
    
    
    print("Writing output file: %s" % (args.output_filename))
    geno_df.to_csv(args.output_filename, sep='\t', index = False, header = True, compression='gzip')


    
    
if __name__ == '__main__':
  main()


