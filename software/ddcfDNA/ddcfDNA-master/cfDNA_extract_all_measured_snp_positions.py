#!/usr/bin/env python3


import sys
import os
import getopt
import argparse
import pandas as pd


def main():
  

  parser = argparse.ArgumentParser("Joining menifest files")
  
  parser.add_argument("-o", "--output_file", dest='output_file', default="", 
                      help=('output_file'))
  
  parser.add_argument("-f", "--manifest_files", dest='manifest_files', default=[], nargs='*',
                      help=('menifest files'))
  
  args = parser.parse_args()
  
  illumina_manifest_files = args.manifest_files
  
  #print(illumina_manifest_files)
  #print(args.output_file)
  
  
  res_df = pd.DataFrame(data=None)
  
  for manifest_file in illumina_manifest_files:
    print("parsing file: %s" % (manifest_file))
    cur_mani_df = pd.read_table(manifest_file, sep='\t')
    cur_mani_df = cur_mani_df[['Chr', 'Position']][cur_mani_df.Chr.isin(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22'])]
    res_df = pd.concat([res_df,cur_mani_df],ignore_index=True)
    res_df.drop_duplicates(inplace=True)

  #res_df.sort_values(inplace = True)
  print("wrting file: %s" % (args.output_file))
  res_df.to_csv(args.output_file, sep='\t', index = False)


if __name__ == '__main__':
  main()