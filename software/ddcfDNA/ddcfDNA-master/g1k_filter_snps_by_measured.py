#!/usr/bin/env python3

"""
Eilon Sharon
"""


from __future__ import division

import sys
import os
import getopt
import argparse
  
import pandas as pd
### SLIGHT MODIFICATION TO CODE: vcf NOT FOUND
import vcf


################################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Filtering 1000 genomes SNPs by meaured SNPs and quality and extracting allele frequency. only single bp REF and ALT are considered")

    # input files

    parser.add_argument("input_snps_vcf_filename", help= "SNPs VCF file")


    parser.add_argument("input_measured_snps_positions_filename", \
                      help="Chr\tPosition tab delimited file")
    
    #output files
    
    parser.add_argument("output_filename", help='tab delimited file of the measured SNPs with G1K information')

    
    
    # options

    parser.add_argument("-p", "--pop", dest='pop', type=str, default='',
                      help='population name, if empty (defualt, will consider the file as all samples AF and super populations. if a population is given it is a filterd file like in http://www.internationalgenome.org/faq/how-can-i-get-allele-frequency-my-variant/)')
    
    parser.add_argument("-q", "--min_quality", dest='min_quality', default=90, type=int,
                      help=('min quality of the SNP in G1K'))
    
    args = parser.parse_args()
    
    ################################################################################################################


    # running the function
    
    print('loading measured SNPs file: %s' % (args.input_measured_snps_positions_filename))
    measured_SNPs_df = pd.read_table(args.input_measured_snps_positions_filename, sep='\t')
    
    print('parsing measured SNPs file ->  dict')
    measured_SNPs_dict = {}

    for idx,row in measured_SNPs_df.iterrows():
        if idx % 200000 == 0:
            print("parsing SNP line %d" %(idx))
        
        cur_chr = row['Chr']
        cur_position = row['Position']
        
        if cur_chr not in measured_SNPs_dict:
             measured_SNPs_dict[cur_chr] = {}
        measured_SNPs_dict[cur_chr][cur_position] = True
    
    ################################################################################################################

    print('parsing VCF file')
    
    vcf_reader = vcf.Reader(filename=args.input_snps_vcf_filename)
    cnt=0
    i = 0
    
    if (args.pop == ''):
      print('parsing as all samples file')
      
      with open(args.output_filename, 'w') as f_out:
          f_out.write("Chr\tPosition\tg1000_ID\tStrand\tg1000_Allele1\tg1000_Allele2\tAF\t" + \
                      "EAS_AF\tAMR_AF\tAFR_AF\tEUR_AF\tSAS_AF\tDP\tAN\tNS\tAC\n")
          for record in vcf_reader:
              i += 1
              if i % 50000 == 0:
                  print("record %d, pass %d " % (i,cnt))
              if ( (len(record.REF) == 1) and
                   (len(record.ALT) == 1) and (len(str(record.ALT[0])) == 1) and
                   (record.QUAL >= args.min_quality) and
                   ( ('chr' + record.CHROM) in measured_SNPs_dict) and
                   (record.POS in measured_SNPs_dict['chr' + record.CHROM])):
                  cnt += 1
                  
                  f_out.write("%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n" % (
                             'chr' + record.CHROM,
                             record.POS, record.ID, '+', record.REF, record.ALT[0],
                             record.INFO['AF'][0],record.INFO['EAS_AF'][0],record.INFO['AMR_AF'][0],
                             record.INFO['AFR_AF'][0],record.INFO['EUR_AF'][0],record.INFO['SAS_AF'][0],
                             record.INFO['DP'], record.INFO['AN'], record.INFO['NS'],record.INFO['AC'][0]))
              
    else:
      print('parsing as single population file: %s' % (args.pop))
      
      with open(args.output_filename, 'w') as f_out:
        f_out.write("Chr\tPosition\tg1000_ID\t%s_AF\t%s_AC\t%s_AN\n" % (args.pop,args.pop,args.pop))
        for record in vcf_reader:
            i += 1
            if i % 50000 == 0:
                print("record %d, pass %d " % (i,cnt))
            if ( (len(record.REF) == 1) and
                 (len(record.ALT) == 1) and (len(str(record.ALT[0])) == 1) and
                 (record.QUAL > args.min_quality) and
                 ( ('chr' + record.CHROM) in measured_SNPs_dict) and
                 (record.POS in measured_SNPs_dict['chr' + record.CHROM])):
                cnt += 1
                
                f_out.write("%s\t%d\t%s\t%f\t%d\t%d\n" % (
                           'chr' + record.CHROM, record.POS, record.ID,
                            ((record.INFO['AC'][0] + 0.0)/record.INFO['AN']),
                            record.INFO['AC'][0], record.INFO['AN']))
      
    print("Summary - parsed %d SNPs, %d passed filter (%f)"  % (i,cnt,cnt/i) )
      

    
    
if __name__ == '__main__':
  main()




