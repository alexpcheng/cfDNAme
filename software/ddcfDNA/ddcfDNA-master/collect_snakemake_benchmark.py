#!/usr/bin/env python3


"""
Eilon Sharon 2017
"""

import argparse
import glob
import os
import pandas as pd
import numpy as np


def collect_benchmark(benchmark_dir, output_seconds_filename, output_time_filename):

    benchmark_files = glob.glob(benchmark_dir + "/*.*.txt")


    benchmark_sec_df = pd.DataFrame(None)
    benchmark_time_df = pd.DataFrame(None)



    for f_i,cur_filename in enumerate(benchmark_files):

        #print("Reading %s" % (cur_filename))
        cur_benchmark_df = pd.read_table(cur_filename, sep='\t')

        cur_benchmark_df = pd.read_table(cur_filename, sep='\t')
        cur_sample, cur_rule, _ = os.path.basename(cur_filename).split('.')

        print("%d. Sample %s, Rule %s, seconds: %f, time: %s" % (f_i,cur_sample, cur_rule, cur_benchmark_df.ix[0,'s'], cur_benchmark_df.ix[0,'h:m:s']))

        if benchmark_sec_df.shape[0] == 0:
            benchmark_sec_df = pd.DataFrame({'Sample' : [cur_sample], cur_rule : [cur_benchmark_df.ix[0,'s'] ]})
            benchmark_time_df = pd.DataFrame({'Sample' : [cur_sample], cur_rule : [cur_benchmark_df.ix[0,'h:m:s'] ]})
        else:
            if cur_rule not in benchmark_sec_df.columns.tolist():
                benchmark_sec_df[cur_rule] = np.nan
                benchmark_time_df[cur_rule] = ''

            if cur_sample in benchmark_sec_df['Sample'].tolist():
                benchmark_sec_df.ix[benchmark_sec_df['Sample']==cur_sample,cur_rule] = cur_benchmark_df.ix[0,'s']
                benchmark_time_df.ix[benchmark_sec_df['Sample']==cur_sample,cur_rule] = cur_benchmark_df.ix[0,'h:m:s']
            else:
                benchmark_sec_df.loc[benchmark_sec_df.shape[0]] = [None] * benchmark_sec_df.shape[1]
                benchmark_sec_df.ix[benchmark_sec_df.shape[0]-1,'Sample'] = cur_sample
                benchmark_sec_df.ix[benchmark_sec_df.shape[0]-1,cur_rule] = cur_benchmark_df.ix[0,'s']

                benchmark_time_df.loc[benchmark_time_df.shape[0]] = [None] * benchmark_time_df.shape[1]
                benchmark_time_df.ix[benchmark_time_df.shape[0]-1,'Sample'] = cur_sample
                benchmark_time_df.ix[benchmark_time_df.shape[0]-1,cur_rule] = cur_benchmark_df.ix[0,'h:m:s']


    print('Saving the benchmark tables to: %s and %s' % (output_seconds_filename, output_time_filename))
    benchmark_sec_df.to_csv(output_seconds_filename, sep='\t', index = False, header = True)
    benchmark_time_df.to_csv(output_time_filename, sep='\t', index = False, header = True)



###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Collecting_snakemake_benchmarks. Files need to have the structure <sample/item>.<rule tag>.txt")

    # input files

    parser.add_argument("benchmark_dir", type=str, help= "snakemake benchmark directory")

    #output files
    parser.add_argument("output_seconds_filename", type=str, 
                      help="output tab-separated values (tsv) file benchmark in seconds")
    
    parser.add_argument("output_time_filename", type=str, 
                      help="output tab-separated values (tsv) file benchmark in hh:mm:ss")
     
    args = parser.parse_args()
    
    collect_benchmark(args.benchmark_dir, args.output_seconds_filename, args.output_time_filename)


if __name__ == '__main__':
  main()


