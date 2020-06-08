#!/usr/bin/env python3

"""
Merge two samples genotype files to create a simulation genotype file
"""

import os
import argparse
import gzip

import itertools as iter

import numpy as np
import pandas as pd



###############################################################################################################
# functions for IBD
###############################################################################################################


def cal_recombination_prob(genetic_distances):
    """
    calculates the probability of an odd number of recombination events
    """
    theta_recomb_p = (1 - np.exp(-2*(genetic_distances/100)))/ 2
    return(theta_recomb_p)


def cal_ibd_y1_vec(genetic_distances,m):
    theta_recomb_p = cal_recombination_prob(genetic_distances)
    return( (1-theta_recomb_p)**(m-2) )

def cal_ibd_y2_vec(genetic_distances):
    theta_recomb_p = cal_recombination_prob(genetic_distances)
    return( (theta_recomb_p)**2 + (1-theta_recomb_p)**2 )


def make_l_sub_transition_matrix(genetic_distance,m,y1,y2):
    """
    transition matrix for one pair of recipient-donor chromosomes
    """
    lsub_trans_mat = np.log(np.array( [[ 1 - ( (1-y1*y2) / (2**(m-1) - 1))  , (1-y1*y2) / (2**(m-1) - 1) ],
                                       [ 1 - y1*y2                          , y1*y2                      ]]) )
    
    # TODO remove
    # in case the assumption that there aere at least two meiosis events does not hold
    # then it can only be twin or parent-child. in both cases there is no recombination and there is always IBD=1
    if (m < 2.0):
        if np.any(np.logical_or(sub_trans_mat<0,sub_trans_mat>1)):
            lsub_trans_mat = np.array( [[-10000,0],
                                       [-10000,0]])
    elif (m > 20): #DEBUG
        lsub_trans_mat = np.array( [[0,-10000],
                                   [0,-10000]])

    return(lsub_trans_mat)

def make_l_transition_matrix(genetic_distance,m_1,m_2,
                           m1_y1,m2_y1,
                           y2,
                           inf_genetic_distance = 1e6):
    
    if np.isfinite(genetic_distance) and  genetic_distance < inf_genetic_distance:
        
        lsub_mat1 = make_l_sub_transition_matrix(genetic_distance,m_1,m1_y1,y2)
        lsub_mat2 = make_l_sub_transition_matrix(genetic_distance,m_2,m2_y1,y2)
        
        # / 2 because intial state does not descriminate between which pair is in IBD in IBD 1
        #trans_mat = np.array(
        #                [[ sub_mat1[0,0]*sub_mat2[0,0], 
        #                   sub_mat1[0,0]*sub_mat2[0,1] + sub_mat1[0,1]*sub_mat2[0,0] , 
        #                   sub_mat1[0,1]*sub_mat2[0,1] ],
        #                 [ (sub_mat1[0,0]*sub_mat2[1,0] + sub_mat1[1,0]*sub_mat2[0,0])/2.0 ,  
        #                   (sub_mat1[0,0]*sub_mat2[1,1] + sub_mat1[0,1]*sub_mat2[1,0] + sub_mat1[1,0]*sub_mat2[0,1] + sub_mat1[1,1]*sub_mat2[0,0])/2.0,
        #                   (sub_mat1[0,1]*sub_mat2[1,1] + sub_mat1[1,1]*sub_mat2[0,1])/2.0],
        #                 [ sub_mat1[1,0]*sub_mat2[1,0], 
        #                   sub_mat1[1,0]*sub_mat2[1,1] + sub_mat1[1,1]*sub_mat2[1,0] , 
        #                   sub_mat1[1,1]*sub_mat2[1,1] ]])
        
        ltrans_mat = np.array(
                        [[ lsub_mat1[0,0]+lsub_mat2[0,0], 
                           lsub_mat1[0,0]+lsub_mat2[0,1],
                           lsub_mat1[0,1]+lsub_mat2[0,0] , 
                           lsub_mat1[0,1]+lsub_mat2[0,1] ],
                         [ lsub_mat1[0,0]+lsub_mat2[1,0],
                           lsub_mat1[0,0]+lsub_mat2[1,1],
                           lsub_mat1[0,1]+lsub_mat2[1,0],
                           lsub_mat1[0,1]+lsub_mat2[1,1] ],
                         [ lsub_mat1[1,0]+lsub_mat2[0,0],
                           lsub_mat1[1,0]+lsub_mat2[0,1],
                           lsub_mat1[1,1]+lsub_mat2[0,0],
                           lsub_mat1[1,1]+lsub_mat2[0,1] ],
                         [ lsub_mat1[1,0]+lsub_mat2[1,0], 
                           lsub_mat1[1,0]+lsub_mat2[1,1],
                           lsub_mat1[1,1]+lsub_mat2[1,0] , 
                           lsub_mat1[1,1]+lsub_mat2[1,1] ] ])
    else:
        
        
        # the use of -10,000 instead of -np.inf is because autograd has probelm dealing with the inf
        
        if (m_1 <= 20):
            marginal_p_noibd_1   = np.log(1-2**(1-m_1))
            marginal_p_ibd_1     = np.log(2**(1-m_1))
        else:
            marginal_p_noibd_1   = 0
            marginal_p_ibd_1     = -10000
        
        if (m_2 <= 20):
            marginal_p_noibd_2   = np.log(1-2**(1-m_2))
            marginal_p_ibd_2     = np.log(2**(1-m_2))
        else:
            marginal_p_noibd_2   = 0
            marginal_p_ibd_2     = -10000
        
        ltrans_mat = np.array(
                        [[marginal_p_noibd_1+marginal_p_noibd_2 , marginal_p_noibd_1+marginal_p_ibd_2 , marginal_p_ibd_1+marginal_p_noibd_2, marginal_p_ibd_1+marginal_p_ibd_2],
                         [marginal_p_noibd_1+marginal_p_noibd_2 , marginal_p_noibd_1+marginal_p_ibd_2 , marginal_p_ibd_1+marginal_p_noibd_2, marginal_p_ibd_1+marginal_p_ibd_2],
                         [marginal_p_noibd_1+marginal_p_noibd_2 , marginal_p_noibd_1+marginal_p_ibd_2 , marginal_p_ibd_1+marginal_p_noibd_2, marginal_p_ibd_1+marginal_p_ibd_2],
                         [marginal_p_noibd_1+marginal_p_noibd_2 , marginal_p_noibd_1+marginal_p_ibd_2 , marginal_p_ibd_1+marginal_p_noibd_2, marginal_p_ibd_1+marginal_p_ibd_2]])
    
    return(ltrans_mat)


def add_ibd_state_per_block(merged_sample_df, ibd_m1, ibd_m2):
    """
    calculating IBD state per block
    """
    
    # pre-cal values
    obs = {}
    obs['inf_genetic_distance'] = 1e6
    obs['genetic_blocks'] = merged_sample_df.genetic_blocks.unique()
    obs['num_blocks'] = obs['genetic_blocks'].shape[0]
    tmp_genetic_distances = merged_sample_df['genetic_avg_distance_FromPrevBlock'].groupby(merged_sample_df['genetic_blocks']).mean().values
    tmp_genetic_distances[np.isinf(tmp_genetic_distances)] = obs['inf_genetic_distance']
    obs['genetic_distances'] = tmp_genetic_distances
    obs['y2_vec']  = cal_ibd_y2_vec(obs['genetic_distances'])
    
    # cal IBD (sampling from distribution)
    merged_sample_df['IBDstate'] = 0

    inf_genetic_distance = obs['inf_genetic_distance']
    genetic_distances = obs['genetic_distances']
    y2_vec = obs['y2_vec']
    m1_y1_vec = cal_ibd_y1_vec(genetic_distances,ibd_m1)
    m2_y1_vec = cal_ibd_y1_vec(genetic_distances,ibd_m2)
    
    prev_state_p = np.array([1.0,1.0,1.0,1.0])/4.0
    
    ibd_state_cnt0 = 0
    ibd_state_cnt1 = 0
    ibd_state_cnt2 = 0
    ibd_state_cnt3 = 0
    
    for i,b in enumerate(obs['genetic_blocks']):
        
        trans_mat = np.exp(make_l_transition_matrix(genetic_distances[i],ibd_m1,ibd_m2,
                                               m1_y1_vec[i],m2_y1_vec[i],
                                               y2_vec[i],
                                               inf_genetic_distance = inf_genetic_distance))

        cur_state_p =  np.sum( np.transpose(np.array([prev_state_p])) * trans_mat,axis=0)
    
        cur_ibd_state = np.sum(np.cumsum(cur_state_p) < np.random.uniform())
        cur_state_p[:]=0.0
        cur_state_p[cur_ibd_state] = 1.0
        prev_state_p = cur_state_p
        
        cur_ibd_state012 = cur_ibd_state
        
        if (cur_ibd_state == 0):
          ibd_state_cnt0 = ibd_state_cnt0 +1
        elif (cur_ibd_state == 1):
          ibd_state_cnt1 = ibd_state_cnt1 +1
        elif (cur_ibd_state == 2):
          ibd_state_cnt2 = ibd_state_cnt2 +1
        elif (cur_ibd_state == 3):
          ibd_state_cnt3 = ibd_state_cnt3 +1
        
        if cur_ibd_state == 1 or cur_ibd_state == 2:
            cur_ibd_state012 = 1
        elif cur_ibd_state == 3:
            cur_ibd_state012 = 2
        
        merged_sample_df.ix[merged_sample_df.genetic_blocks == b, 'IBDstate'] =  cur_ibd_state012
    print("ibd counters:  %d, %d, %d, %d" % (ibd_state_cnt0,ibd_state_cnt1,ibd_state_cnt2,ibd_state_cnt3))
    print("Fraction of SNPs with IBD 0: %f" % (np.sum(merged_sample_df['IBDstate'] == 0.0) / (1.0*merged_sample_df.shape[0])))
    print("Fraction of SNPs with IBD 1: %f" % (np.sum(merged_sample_df['IBDstate'] == 1.0) / (1.0*merged_sample_df.shape[0])))
    print("Fraction of SNPs with IBD 2: %f" % (np.sum(merged_sample_df['IBDstate'] == 2.0) / (1.0*merged_sample_df.shape[0])))

    return(merged_sample_df)

###############################################################################################################
# functions for sim
###############################################################################################################


def select_SNPs_and_reads(merged_sample_df,sample2_fraction, sim_num_snps, sim_snp_selection_method,reads_per_snp):
    """
    selecting SNPs and reads that will be used for the simulation
    """
    
    # selecting SNPs
    merged_sample_df['selected_sim_SNP'] = False
    if sim_snp_selection_method == 'Random':
        sel_snps = np.random.permutation(range(merged_sample_df.shape[0]))[0:min(merged_sample_df.shape[0], sim_num_snps)]
    elif sim_snp_selection_method == 'KMostCommonInEachBlock':
        sel_snps = np.sort(merged_sample_df.groupby(['genetic_blocks'])['G1K_AF'].nlargest(sim_num_snps).index.levels[1])
    else:
        raise ValueError("unknown SNP selection method: |%s|" % (sim_snp_selection_method))
    merged_sample_df.ix[sel_snps,'selected_sim_SNP'] = True
    
    print("Number of selected SNPs: %d" % (merged_sample_df['selected_sim_SNP'].sum()))
    

    # selecting reads for the simulation
    reads_in_sel_snps_sample1 = (merged_sample_df['s1_g1000_Allele1_seq_cnt'][merged_sample_df['selected_sim_SNP']].sum() + 
                             merged_sample_df['s1_g1000_Allele2_seq_cnt'][merged_sample_df['selected_sim_SNP']].sum())

    reads_in_sel_snps_sample2 = (merged_sample_df['s2_g1000_Allele1_seq_cnt'][merged_sample_df['selected_sim_SNP']].sum() + 
                                 merged_sample_df['s2_g1000_Allele2_seq_cnt'][merged_sample_df['selected_sim_SNP']].sum())
    
    print("Number of reads in sample 1: %d" %(reads_in_sel_snps_sample1))
    print("Number of reads in sample 2: %d" %(reads_in_sel_snps_sample2))
    
    if sample2_fraction == 1.0:
        s1_total_constrain = 1e16
    else:
        s1_total_constrain = reads_in_sel_snps_sample1/(1.0-sample2_fraction)
    
    if sample2_fraction == 0.0:
        s2_total_constrain = 1e16
    else:
        s2_total_constrain = reads_in_sel_snps_sample2/sample2_fraction
    
    total_sim_num_reads = min(merged_sample_df['selected_sim_SNP'].sum() * reads_per_snp, 
                              min(s1_total_constrain,s2_total_constrain)) 
                              
    print("total number of reads %d" %(total_sim_num_reads))

    sel_num_reads_s1 = int(np.floor(total_sim_num_reads  * (1.0-sample2_fraction) ))
    sel_num_reads_s2 = int(np.floor(total_sim_num_reads  *      sample2_fraction  ))
    
    sel_reads_s1 = np.sort(np.random.permutation(range(reads_in_sel_snps_sample1))[0:sel_num_reads_s1])
    sel_reads_s2 = np.sort(np.random.permutation(range(reads_in_sel_snps_sample2))[0:sel_num_reads_s2])
    
    i_in_s1 = 0
    i_in_s2 = 0
    
    cnt_s1 = 0
    cnt_s2 = 0
    
    merged_sample_df['g1000_Allele1_seq_cnt'] = 0
    merged_sample_df['g1000_Allele2_seq_cnt'] = 0
    
    
    for r, row in merged_sample_df.iterrows():
        
        if r % 10000 == 0:
            print('selecting reads for row %d (out of %d)' % (r, merged_sample_df.shape[0]))
            if (sel_reads_s1.shape[0] > 0) and (i_in_s1 < sel_reads_s1.shape[0]) :
                print('i_in_s1 %d cnt_s1 %d sel_reads_s1[i_in_s1] %d' % (i_in_s1,cnt_s1, sel_reads_s1[i_in_s1]))
            if (sel_reads_s2.shape[0] > 0) and (i_in_s2 < sel_reads_s2.shape[0]):
                print('i_in_s2 %d cnt_s2 %d sel_reads_s2[i_in_s2] %d' % (i_in_s2,cnt_s2, sel_reads_s2[i_in_s2]))
        
        # skipping SNPs that were not selected
        if row['selected_sim_SNP'] == False:
            continue
        
        # adding reads from s1, allele 1
        for _ in range(row['s1_g1000_Allele1_seq_cnt']):
            if i_in_s1 >= sel_reads_s1.shape[0] or sel_reads_s1.shape[0] == 0:
                break
            elif sel_reads_s1[i_in_s1] == cnt_s1: # selected read
                merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] + 1
                i_in_s1 = i_in_s1 + 1
            elif sel_reads_s1[i_in_s1] < cnt_s1:
                raise RuntimeError('Bug in selecting reads 1: %d %d %d %d' % (r,i_in_s1,cnt_s1,sel_reads_s1[i_in_s1] ))
            cnt_s1 = cnt_s1 + 1
            
        # adding reads from s1, allele 2
        for _ in range(row['s1_g1000_Allele2_seq_cnt']):
            if i_in_s1 >= sel_reads_s1.shape[0] or sel_reads_s1.shape[0] == 0:
                break
            elif sel_reads_s1[i_in_s1] == cnt_s1: # selected read
                merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] + 1
                i_in_s1 = i_in_s1 + 1
            elif sel_reads_s1[i_in_s1] < cnt_s1:
                raise RuntimeError('Bug in selecting reads 2: %d %d %d %d' % (r,i_in_s1,cnt_s1,sel_reads_s1[i_in_s1] ))
            cnt_s1 = cnt_s1 + 1
        
        # adding reads from s2, allele 1
        for _ in range(row['s2_g1000_Allele1_seq_cnt']):
            if i_in_s2 >= sel_reads_s2.shape[0] or sel_reads_s2.shape[0] == 0:
                break
            elif sel_reads_s2[i_in_s2] == cnt_s2: # selected read
                if ( (row['IBDstate'] == 0) or (row['IBDstate'] == 1 and (np.random.uniform() < 0.5)) ):
                    merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] + 1
                else:
                    if (np.random.uniform() < row['g1000_Allele1_receiver_p']):
                        merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] + 1
                    else:
                        merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] + 1
                i_in_s2 = i_in_s2 + 1
            elif sel_reads_s2[i_in_s2] < cnt_s2:
                raise RuntimeError('Bug in selecting reads 3: %d %d %d %d' % (r,i_in_s2,cnt_s2,sel_reads_s2[i_in_s2] ))
            cnt_s2 = cnt_s2 + 1
            
        # adding reads from s2, allele 2
        for _ in range(row['s2_g1000_Allele2_seq_cnt']):
            if i_in_s2 >= sel_reads_s2.shape[0]  or sel_reads_s2.shape[0] == 0:
                break
            elif sel_reads_s2[i_in_s2] == cnt_s2: # selected read
                if ( (row['IBDstate'] == 0) or (row['IBDstate'] == 1 and (np.random.uniform() < 0.5)) ):
                    merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] + 1
                else:
                    if (np.random.uniform() < row['g1000_Allele1_receiver_p']):
                        merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele1_seq_cnt'] + 1
                    else:
                        merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] = merged_sample_df.ix[r, 'g1000_Allele2_seq_cnt'] + 1
                i_in_s2 = i_in_s2 + 1
            elif sel_reads_s2[i_in_s2] < cnt_s2:
                raise RuntimeError('Bug in selecting reads 4: %d %d %d %d' % (r,i_in_s2,cnt_s2,sel_reads_s2[i_in_s2] ))
            cnt_s2 = cnt_s2 + 1

  
    return(merged_sample_df)


def merge_samples(sample1_df,sample2_df,snp_identifying_cols):
    """
    Merging smaples for simulation
    """
    
    sample1_df.drop_duplicates(inplace=True, subset = snp_identifying_cols, keep='first')
    sample2_df.drop_duplicates(inplace=True, subset = snp_identifying_cols, keep='first')
    
    # filter genotyping/sequencing errors where possible
    sample1_df = sample1_df.ix [ ~((sample1_df.g1000_Allele1_receiver_p > 0.9) & (sample1_df.g1000_Allele2_seq_cnt > 0 )) ,:]
    sample1_df = sample1_df.ix [ ~((sample1_df.g1000_Allele2_receiver_p > 0.9) & (sample1_df.g1000_Allele1_seq_cnt > 0 )) ,:]
    
    sample2_df = sample2_df.ix [ ~((sample2_df.g1000_Allele1_receiver_p > 0.9) & (sample2_df.g1000_Allele2_seq_cnt > 0 )) ,:]
    sample2_df = sample2_df.ix [ ~((sample2_df.g1000_Allele2_receiver_p > 0.9) & (sample2_df.g1000_Allele1_seq_cnt > 0 )) ,:]
    
    #renaming columns
    sample1_df = sample1_df.rename(index=str, columns={"g1000_Allele1_seq_cnt": "s1_g1000_Allele1_seq_cnt", "g1000_Allele2_seq_cnt": "s1_g1000_Allele2_seq_cnt"})
    sample2_df = sample2_df.rename(index=str, columns={"g1000_Allele1_seq_cnt": "s2_g1000_Allele1_seq_cnt", "g1000_Allele2_seq_cnt": "s2_g1000_Allele2_seq_cnt"})

    
    sample2_slim_df  = sample2_df[snp_identifying_cols + ['s2_g1000_Allele1_seq_cnt', 's2_g1000_Allele2_seq_cnt']].copy()
    
    merged_sample_df = sample1_df.merge(sample2_slim_df, on = snp_identifying_cols, how='inner')
 
    return(merged_sample_df)
  


###############################################################################################################
# main function
###############################################################################################################

def main():

    __version__ = "0.37"
  
    parser = argparse.ArgumentParser("Merge two samples genotype files to create a simulation genotype file. The file needs to be filttered before it can be used for inference.")


    # input files
    parser.add_argument("input_sample1_filename", type=str, help= "tab-separated values (tsv) file. Sample genotype file")
    parser.add_argument("input_sample2_filename", type=str, help= "tab-separated values (tsv) file. Sample genotype file")


    #output files
    
    parser.add_argument("output_merged_samples_filename", type=str, 
                      help="tab-separated values (tsv) file. Sample genotype file, needs to be filttered before it can be used for inference")
    
    
    # options

    parser.add_argument("-f", "--sample2_fraction", dest='sample2_fraction', default=0.5, type=float,
                      help='fraction of sample 2')
    
    parser.add_argument("-s", "--sim_num_snps", dest='sim_num_snps', default=2000000, type=int,
                      help='number of SNPs (the min between this number and the number of SNPs in samples intersect will be used)')
    
    parser.add_argument("-m", "--sim_snp_selection_method", dest='sim_snp_selection_method', type=str, default='Random',
                      help="""
                          snp selection method (default - Random).
                          Random - select sim_num_snps SNP by random
                          KMostCommonInEachBlock - select sim_num_snps SNPs from each block with highest MAF in G1K_AF 
                          """)
    
    parser.add_argument("-r", "--reads_per_snp", dest='reads_per_snp', default=10.0, type=float,
                      help='number of reads per snp (the min between this number * snps and the number of reads require from each sample will be used)')
    
    parser.add_argument("-m1", "--ibd_m1", dest='ibd_m1', default=100, type=int,
                      help='ibd meiosis of the first pair of chromosomes. int >= 2, default 100 (m1 > 20 is considerd no ibd ). Use only for pure sample (before transplant)!')
    parser.add_argument("-m2", "--ibd_m2", dest='ibd_m2', default=100, type=int,
                      help='ibd meiosis of the first pair of chromosomes. int >= 2, default 100 (m2 > 20 is considerd no ibd ). Use only for pure sample (before transplant)!')

    #parser.add_argument("-k", "--k_snp_to_select", dest='k_snp_to_select', default=100, type=int,
    #                  help='number of SNPs to select (positive int). The use depends on the selection method. For Random it is not used. For KMostCommonInEachBlock it is the number of SNPs from each block')

    parser.add_argument('--version', action='version', version='%(prog)s Version:{version}'.format(version=__version__))
    
    args = parser.parse_args()
    
    print("Args:")
    print(args)
    print('-'*40)
    
    assert(args.sample2_fraction >= 0.0 and args.sample2_fraction <= 1.0)
    assert(args.reads_per_snp >= 0.0)
    assert(args.ibd_m1 >= 2 and args.ibd_m1 >= 2)
    #assert(args.k_snp_to_select > 0)

    # running the function
    print("Running version %s" %(__version__))
    
    snp_identifying_cols = ['Chr', 'Position', 'Strand', 'g1000_Allele1', 'g1000_Allele2']
    
    
    print('Loading file: %s' % (args.input_sample1_filename))
    sample1_df = pd.read_table(args.input_sample1_filename, '\t')
    
    print('Loading file: %s' % (args.input_sample2_filename))
    sample2_df = pd.read_table(args.input_sample2_filename, '\t')
    
    print('Merging sample tables: %s, %s' % (args.input_sample1_filename, args.input_sample2_filename))
    merged_sample_df = merge_samples(sample1_df,sample2_df,snp_identifying_cols)
    
    if (args.ibd_m1 <= 20  or args.ibd_m2 <= 20):
      print('Calculating IBD for each block')
      merged_sample_df = add_ibd_state_per_block(merged_sample_df, args.ibd_m1, args.ibd_m2)
    else:
      print('Setting IBD to zero')
      merged_sample_df['IBDstate'] = 0
    
    print("selecting SNPs and reads for simulation")
    merged_sample_df = select_SNPs_and_reads(merged_sample_df,args.sample2_fraction, args.sim_num_snps, args.sim_snp_selection_method, args.reads_per_snp)
  
    
    print("Writing output file: %s" % (args.output_merged_samples_filename))
    merged_sample_df.to_csv(args.output_merged_samples_filename, sep='\t', index = False, header = True, compression='gzip')
    print('Done!')

    
    
if __name__ == '__main__':
  main()
