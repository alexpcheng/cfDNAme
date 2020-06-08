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

def collect_donor_cfDNA_infenerce_results(samples_df_filename, results_files_dir, output_filename):
    
    print("loading samples table: %s" % (samples_df_filename))
    samples_df = pd.read_table(samples_df_filename, sep='\t')
    
    samples_df['OneGenome_IsRunOK'] = False

    samples_df['OneGenome_Sample'] = ''
    samples_df['OneGenome_Population'] =''
    samples_df['OneGenome_P_donor'] = np.nan
    samples_df['OneGenome_P_recipient'] = np.nan
    samples_df['OneGenome_Sequencing_error'] = np.nan
    samples_df['OneGenome_Genotyping_error'] = np.nan
    samples_df['OneGenome_P_IBD_pair1'] = np.nan
    samples_df['OneGenome_P_IBD_pair2'] = np.nan
    samples_df['OneGenome_nll'] = np.nan
    samples_df['OneGenome_number_of_SNPs'] = np.nan
    samples_df['OneGenome_number_of_reads_mapped_to_SNPs'] = np.nan
    samples_df['OneGenome_number_of_genetic_blocks'] = np.nan
    samples_df['OneGenome_input_relatedness'] = np.nan
    samples_df['OneGenome_num_optimization_starting_points'] = np.nan
    samples_df['OneGenome_minimal_allele_frequency'] = np.nan
    
    print('Loading result files:')

    for sample in samples_df.SampleUniqueStr.tolist():
      cur_sample_filename = results_files_dir + '/' + sample + ".param_fit_results.tsv"
      if os.path.exists(os.path.expanduser(cur_sample_filename)):
          print('Collecting results for sample: %s' % (sample))
          cur_res_df = pd.read_table(cur_sample_filename, sep='\t')
  
          for col in cur_res_df.columns:
              samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_' + col] = \
                          cur_res_df.ix[cur_res_df.IsDonorPop == 1,col].values[0]
          
          samples_df.ix[samples_df.SampleUniqueStr == sample, 'OneGenome_IsRunOK'] = True
  
      else:
          print("Sample %s results file (%s) does not exist" % (sample,cur_sample_filename)) 

    #backward compatibility (TODO - remove)
    samples_df['Sample'] = samples_df['OneGenome_Sample']
    samples_df['OneGenome_DonorPop']= samples_df['OneGenome_Population']
    samples_df['OneGenome_DonorFrac'] = samples_df['OneGenome_P_donor']
    samples_df['OneGenome_nLL'] = samples_df['OneGenome_nll']
    
    samples_df['OneGenome_SeqErr'] = samples_df['OneGenome_Sequencing_error']
    samples_df['OneGenome_GenoErrDependent'] = 0
    samples_df['OneGenome_GenoErrIndependent'] = samples_df['OneGenome_Genotyping_error']
    samples_df['OneGenome_param5'] = samples_df['OneGenome_P_IBD_pair1']
    samples_df['OneGenome_param6'] = samples_df['OneGenome_P_IBD_pair2']
    
    samples_df['OneGenome_P_ibd0'] = (1-samples_df['OneGenome_P_IBD_pair1']) * (1-samples_df['OneGenome_P_IBD_pair2'])
    samples_df['OneGenome_P_ibd2'] = samples_df['OneGenome_P_IBD_pair1'] * samples_df['OneGenome_P_IBD_pair2']
    samples_df['OneGenome_P_ibd1'] = 1-samples_df['OneGenome_P_ibd0']-samples_df['OneGenome_P_ibd2']
    
    samples_df['OneGenome_TotalSNPReadsCnt'] = samples_df['OneGenome_number_of_reads_mapped_to_SNPs']
    samples_df['OneGenome_SeqErrLearnStr'] = 'Independent'
    samples_df['OneGenome_GenoErrLearnStr'] = 'Independent'
    samples_df['OneGenome_RelatednessLearnStr'] = samples_df['OneGenome_input_relatedness']
    
    print('Saving output to: %s' % (output_filename))
    samples_df.to_csv(output_filename, sep='\t', index = False, header = True)


###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Collect donor-cfDNA level inference results")

    # input files

    parser.add_argument("samples_df_filename", type=str, help= "Samples tab-separated values (tsv) file. Must contain a SampleUniqueStr column that contains the sample id")

    parser.add_argument("results_files_dir", type=str, 
                      help="a directory that contains the tab-separated values (tsv) results files")
    
    #output files
    
    parser.add_argument("output_filename", type=str, 
                      help="output tab-separated values (tsv) file name")
    
    
    args = parser.parse_args()
    
    collect_donor_cfDNA_infenerce_results(args.samples_df_filename, args.results_files_dir, args.output_filename)
    
    


    
    
if __name__ == '__main__':
  main()


#!/usr/bin/env python

