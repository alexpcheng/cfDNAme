#!/usr/bin/env python3

"""
Eilon Sharon 2017
"""

import os
import argparse
import gzip

import pandas as pd
import numpy as np

def merge_genotypes_and_alelle_cnts(allele_freq_filename, allele_cnt_filename, output_filename, min_SNP_GC_score):
    
    print('loading input tables')
    allele_freq_df = pd.read_table(allele_freq_filename, sep='\t')
    allele_cnt_df = pd.read_table(allele_cnt_filename, sep='\t')
    
    print("Allele frequencies number of SNPs: %d" % (allele_freq_df.shape[0]))
    print("Allele count       number of SNPs: %d" % (allele_cnt_df.shape[0]))
    
    print('Merging the tables')
    geno_df = pd.merge(allele_freq_df, allele_cnt_df, on = ['Chr','Position'])
    geno_df.drop_duplicates(inplace=True)
    print("Number of SNPs: %d" % (geno_df.shape[0]))
    
    initial_number_of_snps = geno_df.shape[0]
    
    #######################################################################
    # filtering the data
    #######################################################################
    
    print('1. Filtering SNPS with GCScore < %f, passed filter:  %d (filttered out %d)' % (min_SNP_GC_score,(geno_df.GCScore >=min_SNP_GC_score).sum(),(geno_df.GCScore <min_SNP_GC_score).sum() ))
    geno_df = geno_df.ix[geno_df.GCScore >=min_SNP_GC_score,]
    
    print('2. Filtering SNPS with in/dels, passed filter:  %d (filttered out %d)' % ((geno_df.NotSNP == 0).sum(), (geno_df.NotSNP != 0).sum()) )
    geno_df = geno_df.ix[geno_df.NotSNP == 0,]
    
    print('3. Filtering SNPs with no reads, passed filter:  %d (filttered out %d)' % ( (geno_df.SNP_Reads > 0).sum(),(geno_df.SNP_Reads <= 0).sum()  ))
    geno_df = geno_df.ix[geno_df.SNP_Reads > 0,]
    
    ok_geno_filter = ( ( (geno_df.Allele1 == geno_df.g1000_Allele1) | 
                         (geno_df.Allele1 == geno_df.g1000_Allele2) ) &
                       ( (geno_df.Allele2 == geno_df.g1000_Allele1) | 
                         (geno_df.Allele2 == geno_df.g1000_Allele2) ) )
                         
    print('4. Filtering SNPs in whch recipient genotype is not one of the population alleles, %d (filttered out %d)' % ((ok_geno_filter).sum(),(~ok_geno_filter).sum()))
    geno_df = geno_df.ix[ok_geno_filter,]
    


    ok_reads_filter = ( ( (geno_df['A'] == 0) | ("A" == geno_df.g1000_Allele1) | ("A" == geno_df.g1000_Allele2) ) &
                        ( (geno_df['G'] == 0) | ("G" == geno_df.g1000_Allele1) | ("G" == geno_df.g1000_Allele2) ) &
                        ( (geno_df['C'] == 0) | ("C" == geno_df.g1000_Allele1) | ("C" == geno_df.g1000_Allele2) ) &
                        ( (geno_df['T'] == 0) | ("T" == geno_df.g1000_Allele1) | ("T" == geno_df.g1000_Allele2) ) )
    
    
    print('5. Filter SNPs that contain reads that do not map to any of the two alleles, %d (filttered out %d)' % ((ok_reads_filter).sum(),(~ok_reads_filter).sum()))
    
    geno_df = geno_df.ix[ok_reads_filter,]
    
    chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
            'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

    print('6. Filter SNPS that are not in autosomal chromosomes, %d (filttered out %d)' % ( (geno_df.Chr.isin(chrs)).sum(), (~geno_df.Chr.isin(chrs)).sum()))

    geno_df = geno_df.ix[geno_df.Chr.isin(chrs),]
    
    print("Number of SNPs: %d" % (geno_df.shape[0]))
    
    #######################################################################
    # set genotype and allele cfDNA counts
    #######################################################################
    
    print('Setting genotype and allele count columns')
    
    # summing the stats
    geno_df['g1000_Allele1_receiver_p'] = (geno_df.g1000_Allele1 == geno_df.Allele1) * 0.5 + \
                                          (geno_df.g1000_Allele1 == geno_df.Allele2) * 0.5
    
    geno_df['g1000_Allele2_receiver_p'] = (geno_df.g1000_Allele2 == geno_df.Allele1) * 0.5 + \
                                          (geno_df.g1000_Allele2 == geno_df.Allele2) * 0.5
        
    
    geno_df['g1000_Allele1_seq_cnt'] = 0.0
    geno_df['g1000_Allele2_seq_cnt'] = 0.0
    
    geno_df.ix[geno_df.g1000_Allele1 =="A",'g1000_Allele1_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele1 =="A",'A']
    geno_df.ix[geno_df.g1000_Allele1 =="G",'g1000_Allele1_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele1 =="G",'G']
    geno_df.ix[geno_df.g1000_Allele1 =="C",'g1000_Allele1_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele1 =="C",'C']
    geno_df.ix[geno_df.g1000_Allele1 =="T",'g1000_Allele1_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele1 =="T",'T']

    geno_df.ix[geno_df.g1000_Allele2 =="A",'g1000_Allele2_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele2 =="A",'A']
    geno_df.ix[geno_df.g1000_Allele2 =="G",'g1000_Allele2_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele2 =="G",'G']
    geno_df.ix[geno_df.g1000_Allele2 =="C",'g1000_Allele2_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele2 =="C",'C']
    geno_df.ix[geno_df.g1000_Allele2 =="T",'g1000_Allele2_seq_cnt'] = geno_df.ix[geno_df.g1000_Allele2 =="T",'T']

    if ('AF' in geno_df.columns.tolist()):
        geno_df['G1K_AF'] = geno_df['AF']

    AF_cols = [col for col in geno_df.columns.tolist() if "_AF" in col]
    out_cols = ["Chr","Position","g1000_ID","Strand", "GCScore",
                "g1000_Allele1","g1000_Allele2",
                "g1000_Allele1_pop_p","g1000_Allele2_pop_p",
                "g1000_Allele1_receiver_p","g1000_Allele2_receiver_p",
                "g1000_Allele1_seq_cnt","g1000_Allele2_seq_cnt"] + AF_cols
    
    geno_df = geno_df.ix[:,out_cols]
    
    
    #######################################################################
    # printing output table
    #######################################################################
    print('saving the output to: %s' % (output_filename))
    geno_df.to_csv(output_filename, sep='\t', index = False, header = True, compression='gzip')

    num_reads_mapped_to_SNPs = geno_df['g1000_Allele1_seq_cnt'].sum() + geno_df['g1000_Allele2_seq_cnt'].sum()
    print("Number of SNPs: %d" % (geno_df.shape[0]))
    print("Number of Reads: %d" % (num_reads_mapped_to_SNPs))
    print("Avg. coverage of SNP: %f" % (float(num_reads_mapped_to_SNPs)/geno_df.shape[0]))
    print("Fraction of SNPs that pass the filters: %f" % (float(geno_df.shape[0])/initial_number_of_snps) )


###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Mergeing genotyping and population frequencies file with allele counts file and filtering SNPs")

    # input files

    parser.add_argument("allele_freq_filename", type=str, help= "Tab-separated values (tsv) file. Genotyping and populations allele frequencies table")

    parser.add_argument("allele_cnt_filename", type=str, 
                      help="Tab-separated values (tsv) file. Alleles counts (parsed from mpileup results)")
    
    #output files
    
    parser.add_argument("output_filename", type=str, 
                      help="output tab-separated values (tsv) file name")
    
    
    parser.add_argument("-g", "--min_SNP_GC_score", dest='min_SNP_GC_score', default=0.7, type=float,
                      help='minimal SNP genotyping GC score')
    
    
    args = parser.parse_args()
    
    merge_genotypes_and_alelle_cnts(args.allele_freq_filename, args.allele_cnt_filename, args.output_filename, args.min_SNP_GC_score)
        
    
if __name__ == '__main__':
  main()



