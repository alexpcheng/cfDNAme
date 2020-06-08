#!/usr/bin/env python3

"""
Collect donor-cfDNA level inference results
Eilon Sharon 2017
"""

import os
import argparse
import gzip

import pandas as pd
import numpy as np

import pysam

def collect_donor_cfDNA_infenerce_results(samples_df_filename, output_filename, filterred_seq_reads_dir, mapping_dir, wasp_out_dir, sample_genotypes_dir):
    
    print("loading samples table: %s" % (samples_df_filename))
    samples_df = pd.read_table(samples_df_filename, sep='\t')
    
    
    ################################################################
    print('Getting number of Raw reads')
    
    samples_df['OneGenome_NumRawSeqReads'] = np.nan
    for sample in samples_df.SampleUniqueStr.tolist():
        
        print("Sample (raw reads): %s" % (sample))
        
        cur_file_R1 = filterred_seq_reads_dir + '/' + sample + ".R1.filterred.fastq.gz"
        cur_file_R2 = filterred_seq_reads_dir + '/' + sample + ".R2.filterred.fastq.gz"
        
        cur_num_reads = 0.0
        if os.path.exists(os.path.expanduser(cur_file_R1)):
            for line in gzip.open(os.path.expanduser(cur_file_R1),'r'):
                cur_num_reads = cur_num_reads + 1
        if os.path.exists(os.path.expanduser(cur_file_R2)):
            for line in gzip.open(os.path.expanduser(cur_file_R2),'r'):
                cur_num_reads = cur_num_reads + 1
            #cur_num_reads = cur_num_reads + sum(1 for line in gzip.open(os.path.expanduser(cur_file_R2),'r'))
            
        cur_num_reads = cur_num_reads / 4
        samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumRawSeqReads'] = cur_num_reads
    
    ################################################################
    print('Getting number of mapped reads')
    
    samples_df['OneGenome_NumMappedReads'] = np.nan
    
    for sample in samples_df.SampleUniqueStr.tolist():
        print("Sample (mapped reads): %s" % (sample))
        cur_map_bam = mapping_dir + '/' + sample + ".bam"
        
        cur_num_mapped_reads = 0.0
        if os.path.exists(os.path.expanduser(cur_map_bam)):
            samfile = pysam.AlignmentFile("ex1.bam", "rb")
            for x in samfile.fetch():
                cur_num_mapped_reads = cur_num_mapped_reads+1
            
        samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumMappedReads'] = cur_num_reads
    
    
    ################################################################
    print('Getting number of reads mapped to SNPs (after WASP)')
        
    samples_df['OneGenome_NumReadsMapToOptionalSNPs'] = np.nan
    
    for sample in samples_df.SampleUniqueStr.tolist():
        print("Sample (reads mapped to all considered SNPs): %s" % (sample))
        cur_map_bam = wasp_out_dir + '/' + sample + ".bam"
        
        cur_num_mapped_reads = 0.0
        if os.path.exists(os.path.expanduser(cur_map_bam)):
            samfile = pysam.AlignmentFile(cur_map_bam, "rb")
            for x in samfile.fetch():
                cur_num_mapped_reads = cur_num_mapped_reads+1
            
        samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumReadsMapToOptionalSNPs'] = cur_num_reads
    
    
    ################################################################
    print('Getting number of reads mapped to measured SNPs ')
              
    samples_df['OneGenome_NumReadsMapToSNPs'] = np.nan
    samples_df['OneGenome_NumSNPsWithAtLeastOneRead'] = np.nan
    samples_df['OneGenome_NumReadsMapToSNPsAfterFilter'] = np.nan
    samples_df['OneGenome_NumSNPsWithAtLeastOneReadAfterFilter'] = np.nan
    
    
    for sample in samples_df.SampleUniqueStr.tolist():
        print("Sample (reads mapped to measured SNPs): %s" % (sample))
        
        
        cur_genopyte_filename = sample_genotypes_dir + '/' + sample + ".allele_cnt.geno.dist.gz"
        if os.path.exists(os.path.expanduser(cur_genopyte_filename)):
            cur_geno_df = pd.read_table(cur_genopyte_filename, sep='\t')
            samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumReadsMapToSNPs'] = cur_geno_df.shape[0]
            samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumSNPsWithAtLeastOneRead'] = cur_geno_df.g1000_Allele1_seq_cnt.sum() + cur_geno_df.g1000_Allele2_seq_cnt.sum()
        
        cur_filt_genopyte_filename = sample_genotypes_dir + '/' + sample + ".allele_cnt.geno.dist.filt.gz"
        if os.path.exists(os.path.expanduser(cur_filt_genopyte_filename)):
            cur_geno_filt_df = pd.read_table(cur_filt_genopyte_filename, sep='\t')
            samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumReadsMapToSNPsAfterFilter'] = cur_geno_filt_df.shape[0]
            samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_NumSNPsWithAtLeastOneReadAfterFilter'] = cur_geno_filt_df.g1000_Allele1_seq_cnt.sum() + cur_geno_filt_df.g1000_Allele2_seq_cnt.sum()
    
    print('Saving output to: %s' % (output_filename))
    samples_df.to_csv(output_filename, sep='\t', index = False, header = True)


###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Collect donor-cfDNA inference pipeline stats")

    # input files

    parser.add_argument("samples_df_filename", type=str, help= " Samples tab-separated values (tsv) file - either the input file or the file with the infered dd-cfDNA. Must contain a SampleUniqueStr column that contains the sample id")

    #output files
    
    parser.add_argument("output_filename", type=str, 
                      help="output tab-separated values (tsv) file name")
    

    #
    
    parser.add_argument("-f", "--filterred_seq_reads_dir", dest='filterred_seq_reads_dir', default="Output/sequencing", type=str,
                      help='filttered reads directory')
    
    parser.add_argument("-m", "--mapping_dir", dest='mapping_dir', default="Output/mapping", type=str,
                      help='mapped reads directory')
    
    parser.add_argument("-w", "--wasp_out_dir", dest='wasp_out_dir', default="Output/wasp_out", type=str,
                      help='wasp out mapping reads directory')
    
    parser.add_argument("-g", "--sample_genotypes_dir", dest='sample_genotypes_dir', default="Output/SampleGenotypes", type=str,
                      help='sample genotypes directory')
    
    
    args = parser.parse_args()
    
    collect_donor_cfDNA_infenerce_results(args.samples_df_filename, args.output_filename, 
    																			args.filterred_seq_reads_dir, 
    																			args.mapping_dir, 
    																			args.wasp_out_dir, 
    																			args.sample_genotypes_dir)
    
  
if __name__ == '__main__':
  main()


