#!/usr/bin/env python3

"""
Comparison of two genomes and one genome performances
Eilon Sharon 2017
"""

import os
import argparse
import gzip

import itertools

import pandas as pd
import numpy as np

import matplotlib
#matplotlib.style.use('ggplot')
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.optimize import minimize



from sklearn.metrics import roc_auc_score

import seaborn as sns
sns.set(rc={"figure.figsize": (6, 6)})

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()


###############################################################################################################
# R functions for ROC analysis
###############################################################################################################

robjects.r("library(pROC)")
robjects.r("library(ggplot2)")

robjects.r(
"""
roc_compare <- function(response, pred_g1, pred_g2) {
  
  require("pROC")
  
  proc_pl_df = data.frame("response" = response, "pred_g1" = pred_g1, "pred_g2" = pred_g2)
  
  proc_g1 = roc(as.numeric(proc_pl_df$response), as.numeric(proc_pl_df$pred_g1) )
  proc_g2 = roc(as.numeric(proc_pl_df$response), as.numeric(proc_pl_df$pred_g2) )
  proc_test_res = roc.test(proc_g1,proc_g2)
  
  return(proc_test_res$p.value)
  
}
""")


robjects.r(
"""
roc_ci <- function(response, pred) {
  
  require("pROC")
  
  ci = ci.auc(roc(as.numeric(response), as.numeric(pred)), conf.level=0.95)
  return(ci)
  
}
""")


robjects.r(
"""
plot_auc <- function(auc_pl_df, out_filename) {

    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
    limits =  aes(ymax = AUC_down, ymin=AUC_up)
    dodge <- position_dodge(width=0.9)
    
    pdf(file = out_filename[1], useDingbats=FALSE,width=5,height=3)
    pl <- ggplot(auc_pl_df, aes(x=factor(Comparison), y=AUC, fill=factor(Method))) +
    geom_bar(stat="identity",position=dodge) + coord_flip() +
    geom_errorbar(limits, position=dodge, width=0.25) +
    geom_text(aes(label = sapply(as.numeric(auc_pl_df$AUC), function(x) {return( sprintf("%.2f", x) ) }), y = 0.1), size = 3, position=dodge) +
    geom_text(aes(label = sapply(as.numeric(auc_pl_df$pvalue), function(x) {return( sprintf("p-value = %.2f", x) ) }), y = 0.4), size = 3) +
    xlab("Organ rejection") +
    scale_color_manual(values=cbPalette) +
    scale_fill_manual(values=cbPalette) +
    theme_bw(base_size = 10,base_family = "Helvetica") +
    theme(axis.title.x = element_text(size=10,color="black"),
          axis.title.y = element_text(size=10,color="black"),
          axis.text.x=element_text(size=8,color="black") ,
          axis.text.y=element_text(size=8,color="black"),
          legend.position="none",
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(pl)
    dev.off()
}
""")

r_roc_compare = robjects.r['roc_compare']
r_roc_ci= robjects.r['roc_ci']
r_plot_auc= robjects.r['plot_auc']

###############################################################################################################
# R related helper function
###############################################################################################################

def cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, comparison, ind, use_correctedG2 = False):
    roc_pl_df.ix[ind:(ind+1), "Comparison"] = comparison

    roc_ci_g1 = r_roc_ci(robjects.IntVector(response), robjects.FloatVector(pred_g1))
    
    if (use_correctedG2):
        roc_test = r_roc_compare(robjects.IntVector(response), robjects.FloatVector(pred_g1), robjects.FloatVector(pred_g2cor))
        roc_ci_g2 = r_roc_ci(robjects.IntVector(response), robjects.FloatVector(pred_g2cor))
    else:
        roc_test = r_roc_compare(robjects.IntVector(response), robjects.FloatVector(pred_g1), robjects.FloatVector(pred_g2))
        roc_ci_g2 = r_roc_ci(robjects.IntVector(response), robjects.FloatVector(pred_g2))

    roc_pl_df.ix[ind:(ind+1), "pvalue"] = roc_test[0]
    roc_pl_df.ix[ind, "AUC_down"] = roc_ci_g2[0]
    roc_pl_df.ix[ind, "AUC"] = roc_ci_g2[1]
    roc_pl_df.ix[ind, "AUC_up"] = roc_ci_g2[2]
    roc_pl_df.ix[ind+1, "AUC_down"] = roc_ci_g1[0]
    roc_pl_df.ix[ind+1, "AUC"] = roc_ci_g1[1]
    roc_pl_df.ix[ind+1, "AUC_up"] = roc_ci_g1[2]
    
    return(roc_pl_df)


###############################################################################################################
# constants
###############################################################################################################

# TODO
corrupted_sampels = ["L28_W2","L28_M3","L68_M3n1","L68_W5","I5_W10"]


cbPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# The palette with grey:
cbgPalette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# The palette with black:
cbbPalette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# The palette with grey:
cbgPalette = ["#999999", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9"]
cbbPalette = ["#000000", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9"]





###############################################################################################################
# more helper function
###############################################################################################################

def load_data(output_dir, fit_versions, corrupted_sampels, min_2G_seq_error, min_DaysPostTransplant, min_reads_that_map_to_snps):
    """
    loading the data and initial parsing
    """
    res_dfs = {}
    
    for fit_ver in fit_versions:
        cur_res_df = pd.read_table(output_dir + '/Fit/fit_' + fit_ver + '.res.tsv', sep='\t', na_values = "")
        
        print("%s number of samples: %d" %(fit_ver, cur_res_df.shape[0]))
    
        
        # remove corrupted
        cur_res_df = cur_res_df.ix[~cur_res_df.SampleUniqueStr.isin(corrupted_sampels),]
        # remove duplicated lines
        cur_res_df = cur_res_df.ix[~(cur_res_df.SampleUniqueStr.isin(["I3_W3","I3_W4"]) & cur_res_df['Read1'].str.contains("Aligned_reads_131213")) ,]
    
        print("%s number of samples after removing corrupted or duplicated: %d" %(fit_ver, cur_res_df.shape[0]))
    
        
        # filterred samples that were filterred out in the 2G experiment
        cur_res_df = cur_res_df.ix[cur_res_df.Error <= min_2G_seq_error,]
        print("%s number of samples after removing 2G filttered (sequencing error): %d" %(fit_ver, cur_res_df.shape[0]))
        
        cur_res_df = cur_res_df.ix[(cur_res_df.DaysPostTransplant >= min_DaysPostTransplant) & (~cur_res_df.DaysPostTransplant.isnull()),]
        print("%s number of samples after removing 2G filttered (days from transplant): %d" %(fit_ver, cur_res_df.shape[0]))
        
        #print("%s number of samples after removing 2G filttered: %d" %(fit_ver, cur_res_df.shape[0]))
    
        
        cur_res_df = cur_res_df.ix[cur_res_df.OneGenome_TotalSNPReadsCnt >= min_reads_that_map_to_snps,]
        print("%s number of samples after removing sample with small number of reads mapped to SNPs (days from transplant): %d" %(fit_ver, cur_res_df.shape[0]))
        
        
        #print(cur_res_df)
        #print(sum( (~cur_res_df.DaysPostTransplant.isnull())  ))
        #print(sum( (~cur_res_df.DonorCorrected.isnull())  ))
        #print(sum( (~cur_res_df.OneGenome_DonorFrac.isnull())   ))
        #print(sum( (cur_res_df.IsFullInfo == True)  ))
        #print(sum( (cur_res_df.OneGenome_IsRunOK == True)  ))
        
        # filtering partial data
        cur_res_df = cur_res_df.ix[ ( (~cur_res_df.DaysPostTransplant.isnull()) & 
                                      (~cur_res_df.DonorCorrected.isnull()) &
                                      (~cur_res_df.OneGenome_DonorFrac.isnull()) &
                                      (cur_res_df.IsFullInfo == True) &
                                      (cur_res_df.OneGenome_IsRunOK == True) ),]
        
        print("%s number of samples after filtering partial data: %d" %(fit_ver, cur_res_df.shape[0]))
    
        res_dfs[fit_ver] = cur_res_df
        
        if (cur_res_df.shape[0] < 3):
            raise ValueError("Less than 3 sample left after filtering: %s" % (output_dir + '/Fit/fit_' + fit_ver + '.res.tsv'))
        
        
        
    return(res_dfs)

def parse_clinical_terms(fit_versions, res_dfs, set_type):
    """
    parse clinical terms
    """

    Biopsy_order = ["No data"]
    
    for fit_ver in fit_versions:
        out_df = res_dfs[fit_ver]
        out_df['Biopsy_level'] = "No data"
        
        if set_type == "lt":
            
            out_df.Biopsy[out_df.Biopsy.isnull()] = ""
            # parse biopsy
            out_df.ix[out_df.Biopsy.str.startswith('A0'),'Biopsy_level'] = "1. Normal"
            out_df.ix[out_df.Biopsy.str.startswith('A1'),'Biopsy_level'] = "2. Minimal"
            out_df.ix[out_df.Biopsy.str.startswith('A2'),'Biopsy_level'] = "3. Mild"
            out_df.ix[out_df.Biopsy.str.startswith('A3'),'Biopsy_level'] = "4. Moderate"
            out_df.ix[out_df.Biopsy.str.startswith('A4'),'Biopsy_level'] = "5. Severe"
            
            Biopsy_order = ["No data", "1. Normal", "2. Minimal", "3. Mild", "4. Moderate", "5. Severe"]
            
            # parse AMR
            out_df.Clinical_signs_of_AMR[out_df.Clinical_signs_of_AMR.isnull()] = ""
            out_df['AMR'] =  np.nan
            
            out_df.ix[out_df.Clinical_signs_of_AMR.isin(["No","no","NO"]),'AMR'] = 0.0
            out_df.ix[out_df.Clinical_signs_of_AMR.isin(["Yes","yes","YES"]),'AMR'] = 1.0
            out_df.ix[out_df.Clinical_signs_of_AMR.isin(["<NA>", ""]) | out_df.Clinical_signs_of_AMR.isnull(),'AMR'] = np.nan
            
            # parse number of lungs
            out_df['Num_of_Lungs'] = out_df['NofLungs']
    
            # parse treatment
            out_df.Treated_for_rejection[out_df.Treated_for_rejection.isnull()] = ""
            out_df['GotRejectionTreatment'] = np.nan
            out_df.ix[out_df.Treated_for_rejection.isin(["No","no","NO"]),'GotRejectionTreatment'] = 0.0
            out_df.ix[out_df.Treated_for_rejection.isin(["Yes","yes","YES"]),'GotRejectionTreatment'] = 1 
            out_df.ix[out_df.Treated_for_rejection.isin(["<NA>", ""]) | out_df.Treated_for_rejection.isnull(),'GotRejectionTreatment'] = np.nan
        elif set_type == "ht":
            # parse biopsy
            out_df.ix[out_df.Biopsy.isin(["0",0]),'Biopsy_level'] = "1. Quiescence"
            out_df.ix[out_df.Biopsy.isin(["1R/1A","1R/1B","1R/2"]),'Biopsy_level'] = "2. Mild rejection"
            out_df.ix[out_df.Biopsy.isin(["2R/3A","3R/3B","AMR","AMR2"]),'Biopsy_level'] = "3. Moderate or severe rejection"
    
            Biopsy_order = ["No data", "1. Quiescence", "2. Mild rejection", "3. Moderate or severe rejection"]
            
            
            out_df['Is_severe_rejection'] = False
            out_df.ix[out_df.Biopsy.isin(["3R/3B"]),'Is_severe_rejection'] = True
        elif set_type == "bm":
            out_df.Biopsy_level = ""
        else:
            raise ValueError("unknown set name:" + set_type)
        
        res_dfs[fit_ver] = out_df
    return(res_dfs)

def fix_2g_by_seq_error(fit_versions, res_dfs, set_type, fix_by_BC_flag):
    """
    Making sure donor is corrected correctly
    """
    #eps = 1e-3

    for fit_ver in fit_versions:
        res_df = res_dfs[fit_ver]
        
        
        
        
        res_df['DonorCorrected'] = res_df.Donor #- 3.6*res_df.Error
        res_df['DonorCorrectedByError'] = res_df.Donor - 3.6*res_df.Error
        
        #print("%s MAE before: %f" %(fit_ver, np.sum(np.absolute(  (res_df.DonorCorrected) - (100*res_df.OneGenome_DonorFrac)   ))))
        
        #print((res_df.DonorCorrected) - (100*res_df.OneGenome_DonorFrac))
        
        #opt = minimize(lambda x: np.sum(np.absolute(  (res_df.DonorCorrected+x) - (100*res_df.OneGenome_DonorFrac)   )), x0=[0])
        #x = 0 # max(opt.x[0], 1e-5)
        #print("%s shifting G2 by %f (min G1 %f, opt x: %f)" % (fit_ver, x, min((100*res_df.OneGenome_DonorFrac)),opt.x[0]))
        #res_df['DonorFraction_G2'] = x + res_df['DonorCorrected']
        res_df['DonorFraction_G2'] = res_df['DonorCorrected']
        res_df['DonorFraction_G2_CorrectedByError'] = res_df['DonorCorrectedByError']
        
        #print("%s MAE after: %f" %(fit_ver, np.sum(np.absolute(  (res_df.DonorCorrected+opt.x[0]) - (100*res_df.OneGenome_DonorFrac)   ))))
        
        #res_df['DonorFraction_G2_log10'] = np.log10(res_df['DonorFraction_G2']+eps+np.absolute(min(0,res_df['DonorFraction_G2'].min())))
        #res_df['DonorFraction_G1_log10'] = np.log10(res_df['DonorFraction_G1'])
        
        
        if (fix_by_BC_flag and set_type == "bm"):
            
            print('fixing by BC')
                        
            res_bc_df = res_df.ix[res_df.SampleUniqueStr.str.contains('_BC'),['SampleUniqueStr','Patient','OneGenome_DonorFrac']].copy()
            res_bc_df = res_bc_df.rename(index=str, columns={"SampleUniqueStr": "SampleUniqueStr_BC", 
                                                             "OneGenome_DonorFrac" : "OneGenome_DonorFrac_BC"})
            
            res_df = res_df.merge(res_bc_df, how='left', on = "Patient")

            res_df["OneGenome_DonorFrac"] = np.maximum(0.0, res_df["OneGenome_DonorFrac"] - ( (1.0 - res_df["OneGenome_DonorFrac"]) * res_df["OneGenome_DonorFrac_BC"] ).values)
            
            
        res_df['DonorFraction_G1'] = 100*res_df.OneGenome_DonorFrac
        
        
        res_dfs[fit_ver] = res_df
        
    return(res_dfs)


def plot_2G_vs_1G(fit_versions, res_dfs, fig_output_dir, set_type):
    """
    2D plot of 2G vs. 1G
    """
    for fit_ver in fit_versions:
        res_df = res_dfs[fit_ver]
        
        biopsy = res_df.Biopsy_level.unique()
        biopsy.sort()
        biopsy =  list([biopsy[-1]]) + list(biopsy[:-1])
        max_val = max(res_df.DonorFraction_G2.max(), res_df.DonorFraction_G1.max())
        spearman_rho, _ = spearmanr(res_df.DonorFraction_G2, res_df.DonorFraction_G1)
        
        print("MAE: %f" % ( np.mean(np.absolute(res_df.DonorFraction_G2-res_df.DonorFraction_G1)) ))
        
        with sns.axes_style("white"):
            sns.set_style("ticks")
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val*1.1], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val*0.9], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val*1.5], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val/1.5], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val*2], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #ax.plot([0,max_val],[0,max_val*0.5], linestyle=':', color = "gray", label=None,linewidth=0.5)
            #sns.kdeplot(res_df.DonorFraction_G2, res_df.DonorFraction_G1, ax=ax, legend=False, kernel='gau')
            handles = [None] * len(biopsy)
            for b,c in zip(biopsy,cbgPalette[0:len(biopsy)]):
                ax.plot(res_df.DonorFraction_G2[res_df.Biopsy_level == b], res_df.DonorFraction_G1[res_df.Biopsy_level ==b], label=b, marker='.', linestyle='', color=c, markersize=4)
            plt.xscale('log', nonposy='clip')
            plt.yscale('log', nonposy='clip')
            plt.xlabel("Percent of dd-cfDNA estimated using two genomes (log scale)",fontsize=10)
            plt.ylabel("Percent of dd-cfDNA estimated using one genome (log scale)",fontsize=10)
            if (set_type not in ['bm']):
                rejection_legend = plt.legend(loc='upper left')
                rejection_legend.get_title().set_fontsize('10')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            #plt.title(r'Spearman \rho  = %.2f' % (spearman_rho))
            ax.text(0.9, 0.02, r"Spearman's $\rho$  = %.2f" % (spearman_rho),
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax.transAxes,
                    color='black', fontsize=10)
            
            fig.savefig( fig_output_dir + "/" + fit_ver + "_dotplot_log"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None,linewidth=0.5)
            handles = [None] * len(biopsy)
            for b,c in zip(biopsy,cbbPalette[0:len(biopsy)]):
                ax.plot(res_df.DonorFraction_G2[res_df.Biopsy_level == b], res_df.DonorFraction_G1[res_df.Biopsy_level ==b], label=b, marker='.', linestyle='', color=c, markersize=4)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Percent of dd-cfDNA estimated using two genomes",fontsize=10)
            plt.ylabel("Percent of dd-cfDNA estimated using one genome",fontsize=10)
            if (set_type not in ['bm']):
                rejection_legend = plt.legend(loc='upper left',fontsize=8, title="Organ rejection")
                rejection_legend.get_title().set_fontsize('10')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            #plt.title(r'Spearman \rho  = %.2f' % (spearman_rho))
            ax.text(0.9, 0.02, r"Spearman's $\rho$  = %.2f" % (spearman_rho),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes,
            color='black', fontsize=10)
            
            fig.savefig( fig_output_dir + "/" + fit_ver + "_dotplot_lin"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            
            # cumulative absolute error
            abs_error = np.sort(np.absolute(res_df.DonorFraction_G2 - res_df.DonorFraction_G1))
            
            abs_error_cdf = np.array(range(abs_error.shape[0])) * (1.0 / abs_error.shape[0])
            
            re_abs_error = np.sort(abs_error / res_df.DonorFraction_G2)
          
            
              
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot(abs_error, abs_error_cdf, marker='', linestyle='-', color='black', linewidth=1)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Absolute difference between precent of dd-cfDNA\nestimation using one and two genomes methods (\%)",fontsize=10)
            plt.ylabel("Samples with smaller or equal absolute error",fontsize=10)
            #plt.xscale('log', nonposy='clip')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            fig.savefig( fig_output_dir + "/" + fit_ver + "_2g_1G_abs_dif_cdf"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot(re_abs_error, abs_error_cdf, marker='', linestyle='-', color='black', linewidth=1)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Absolute relative difference between precent of dd-cfDNA\nestimation using one and two genomes methods",fontsize=10)
            plt.ylabel("Samples with smaller or equal absolute error",fontsize=10)
            #plt.xscale('log', nonposy='clip')
            plt.xlim((0,2))
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            fig.savefig( fig_output_dir + "/" + fit_ver + "_2g_1G_rel_abs_dif_cdf"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            pearson_R,_ = pearsonr(res_df.DonorFraction_G2, res_df.DonorFraction_G1)
            print('%s spearman rho: %f' % (fit_ver, spearman_rho)) 
            print('%s pearson R^2: %f' % (fit_ver, pearson_R**2)) 
            print('%s median absolute error: %f' % (fit_ver, np.median(abs_error)))
            print('%s median relative absolute error: %f' % (fit_ver, np.median(re_abs_error)))
            
            
            
            
            
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot(res_df.DonorFraction_G2, np.abs((res_df.DonorFraction_G2-res_df.DonorFraction_G1)/res_df.DonorFraction_G2),  marker='.', linestyle='', color='black', markersize = 4)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Percent of dd-cfDNA estimated using two genomes",fontsize=10)
            plt.ylabel("Absolute relative difference between precent of dd-cfDNA\nestimation using one and two genomes methods",fontsize=10)
            plt.xscale('log', nonposy='clip')
            plt.yscale('log', nonposy='clip')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            fig.savefig( fig_output_dir + "/" + fit_ver + "_2g_to_rel_abs_dif_dotplot"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            fig, ax = plt.subplots(1,1,figsize=(4,4)) 
            ax.plot(res_df.DonorFraction_G2, np.abs(res_df.DonorFraction_G2-res_df.DonorFraction_G1),  marker='.', linestyle='', color='black', markersize = 4)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Percent of dd-cfDNA estimated using two genomes",fontsize=10)
            plt.ylabel("Absolute difference between precent of dd-cfDNA\nestimation using one and two genomes methods",fontsize=10)
            plt.xscale('log', nonposy='clip')
            plt.yscale('log', nonposy='clip')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            fig.savefig( fig_output_dir + "/" + fit_ver + "_2g_to_abs_dif_dotplot"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            
            max_val = res_df.DonorFraction_G2.max()
            min_val = 0.01#min(0.01,res_df.DonorFraction_G2.min())
            min_max_vec = np.array([min_val, max_val])
            
            fig, ax = plt.subplots(1,1,figsize=(4,4))
            tf_inds = (res_df.DonorFraction_G2 > 0) & (res_df.DonorFraction_G1 > 0)
            ax.plot(np.log10(min_max_vec),np.log10(0.01*min_max_vec), linestyle=':', color = "gray", label=None, linewidth = 0.5)
            ax.plot(np.log10(min_max_vec),np.log10(0.05*min_max_vec), linestyle=':', color = "gray", label=None, linewidth = 0.5)
            ax.plot(np.log10(min_max_vec),np.log10(0.1*min_max_vec), linestyle=':', color = "gray", label=None, linewidth = 0.5)
            ax.plot(np.log10(min_max_vec),np.log10(0.3*min_max_vec), linestyle=':', color = "gray", label=None, linewidth = 0.5)
            ax.plot(np.log10(min_max_vec),np.log10(1.0*min_max_vec), linestyle=':', color = "gray", label=None, linewidth = 0.5)
            ax.plot(np.log10(res_df.DonorFraction_G2), np.log10(np.abs(res_df.DonorFraction_G2-res_df.DonorFraction_G1)),  marker='.', linestyle='', color='black', markersize = 4)
            ax = sns.kdeplot(np.log10(res_df.DonorFraction_G2[tf_inds]), np.log10(np.abs(res_df.DonorFraction_G2[tf_inds]-res_df.DonorFraction_G1[tf_inds])), ax=ax, shade=True,shade_lowest=False)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            plt.xlabel("Percent of dd-cfDNA estimated using two genomes",fontsize=10)
            plt.ylabel("Absolute difference between percent of dd-cfDNA\nestimation using one and two genomes methods",fontsize=10)
            #plt.xscale('log', nonposy='clip')
            #plt.yscale('log', nonposy='clip')
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter(r"$10^{%d}$")) #ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.FormatStrFormatter(r"$10^{%d}$")) #FuncFormatter(lambda x, pos: r"$10^" + str(x) + "$"))  #FormatStrFormatter("10^%d"))
            plt.xlim(plt.xlim()[0], min(plt.xlim()[1],np.log10(100)))
            fig.savefig( fig_output_dir + "/" + fit_ver + "_2g_to_abs_dif_dotplot_withKDe"  ".pdf" , bbox_inches='tight')
            plt.close()
            
            



def plot_2G_vs_1G_algorithms(fit_versions, res_dfs, fig_output_dir, set_type, models_desc_dict):
    """
    2D plot of 2G vs. 1G algorithms to compare the algorithms
    """
    
    markers_dict = {'o': 'circle', 'v': 'triangle_down', '^': 'triangle_up', '<': 'triangle_left', '>': 'triangle_right', '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right', '8': 'octagon', 's': 'square', 'p': 'pentagon', '*': 'star', 'h': 'hexagon1', 'H': 'hexagon2', '+': 'plus', 'x': 'x', 'D': 'diamond', 'd': 'thin_diamond', '|': 'vline', '_': 'hline', 'P': 'plus_filled', 'X': 'x_filled', 0: 'tickleft', 1: 'tickright', 2: 'tickup', 3: 'tickdown', 4: 'caretleft', 5: 'caretright', 6: 'caretup', 7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase', 10: 'caretupbase', 11: 'caretdownbase'}
    
    markers_vec = list(markers_dict.keys())
    markers_vec = ["o","s","d","+","x","8","h",">","v","<","^","h" ] 
    
    colors = itertools.cycle(cbPalette)
    markers = itertools.cycle(markers_vec)
    #print(markers_vec)
    
    max_val = 0
    for fit_ver in fit_versions:
        res_df = res_dfs[fit_ver]
        max_val = max(max_val,max(res_df.DonorFraction_G2.max(), res_df.DonorFraction_G1.max()))
    
    with sns.axes_style("white"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(4,4))
        ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None)
        ax.plot([0,max_val],[0,max_val/2], linestyle=':', color = "gray", label=None)
        
        res_df = res_dfs[fit_versions[0]]

        # plotting just for patients legend
        markers = itertools.cycle(markers_vec)
        patient_plot_handles = []
        
        for patient in res_df.Patient.unique():
            
            cur_marker = next(markers)
                
            cur_pl, = ax.plot(res_df.DonorFraction_G2[res_df.Patient == patient], res_df.DonorFraction_G1[res_df.Patient == patient], label=str(patient),
                              marker=cur_marker, linestyle='', color='black', fillstyle='none', markeredgecolor = 'black',markeredgewidth=1.0, markersize=4)
            
            patient_plot_handles.append(cur_pl)


        fit_plots_handles = []
        for fit_ver in fit_versions:
            res_df = res_dfs[fit_ver]
            
            cur_color = next(colors)
            markers = itertools.cycle(markers_vec)
            
            patients = res_df.Patient.unique()
            patients.sort()
            
            for p,patient in enumerate(patients):
                
                cur_marker = next(markers)
                    
                cur_pl, = ax.plot(res_df.DonorFraction_G2[res_df.Patient == patient], res_df.DonorFraction_G1[res_df.Patient == patient], label=fit_ver,
                                  marker=cur_marker, linestyle='', color=cur_color, fillstyle='none', markeredgecolor = cur_color,markeredgewidth=1.0, markersize=4)
                
                if (p==0):
                    fit_plots_handles.append(cur_pl)
            

        plt.xlabel("Percent of dd-cfDNA estimated using two genomes",fontsize=10)
        plt.ylabel("Percent of dd-cfDNA estimated using one genome",fontsize=10)
        #plt.legend(loc='upper left') 
        
        legend_patients = plt.legend(patient_plot_handles, res_df.Patient.unique(), loc='upper left',title='Patient',fontsize=8)
        legend_alg = plt.legend(fit_plots_handles, [models_desc_dict[f] for f in fit_versions], loc='lower right', title="Model",fontsize=8)
        
        legend_patients.get_title().set_fontsize('10')
        legend_alg.get_title().set_fontsize('10')

        plt.gca().add_artist(legend_patients)
        plt.gca().add_artist(legend_alg)
        
        sns.despine()
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        
        fig.savefig( fig_output_dir + "/" + fit_ver + "_dotplot_multi_lin"  ".pdf" , bbox_inches='tight')
        plt.close()
            

def plot_patients(fit_versions, res_dfs, fig_output_dir, set_type, models_desc_dict, use_axis_labels = False):
    
    """
    Plot individual patients plots
    """

    marker_fill_styles = itertools.cycle(['left','right','bottom','top'])

    res_df = res_dfs[fit_versions[0]]
    patients = res_df.Patient.unique()
    patients.sort()
    
    for p,patient in enumerate(patients):
        
        print(patient)
        colors = itertools.cycle(cbPalette)
                
        with sns.axes_style("white"):
            sns.set_style("ticks")
            if (use_axis_labels):
                fig, ax = plt.subplots(1,1,figsize=(4,4))
            else:
                fig, ax = plt.subplots(1,1,figsize=(2.5,2.5))
            fit_plots_handles = []
            for f,fit_ver in enumerate(fit_versions):
                res_df = res_dfs[fit_ver]
                
                
                cur_marker = 'o'
                
                cur_res_df = res_df.ix[res_df.Patient == patient, ['DaysPostTransplant', 'DonorFraction_G1','DonorFraction_G2']]
                cur_res_df.sort_values(by=['DaysPostTransplant'], inplace=True)
                
                if (f==0):
                    cur_color = next(colors)
                    cur_marker_fill_style = next(marker_fill_styles)
                    cur_pl, = ax.plot(cur_res_df.DaysPostTransplant, cur_res_df.DonorFraction_G2, label=fit_ver, marker=cur_marker, linestyle='-', color=cur_color, fillstyle=cur_marker_fill_style) # , markeredgecolor = cur_color,markeredgewidth=1.0
                    fit_plots_handles.append(cur_pl)

                cur_color = next(colors)
                cur_marker_fill_style = next(marker_fill_styles)
                cur_pl, = ax.plot(cur_res_df.DaysPostTransplant, cur_res_df.DonorFraction_G1, label=fit_ver, marker=cur_marker, linestyle='-', color=cur_color, fillstyle=cur_marker_fill_style) # , markeredgecolor = cur_color,markeredgewidth=1.0
                fit_plots_handles.append(cur_pl)
            
            if (use_axis_labels):
                plt.xlabel("Days post transplant",fontsize=10)
                plt.ylabel("Estimated percent of dd-cfDNA",fontsize=10)
                
            
            legend_alg = plt.legend(fit_plots_handles, ['Two genomes'] + [models_desc_dict[f] for f in fit_versions], loc='best', title="Model", fontsize=8) #
            legend_alg.get_title().set_fontsize('10')
            
            if (use_axis_labels):
                plt.gca().add_artist(legend_alg)
            
            sns.despine()
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            
            if (use_axis_labels):
                fig.savefig( fig_output_dir + "/Patients/" +  patient +  "_labeledaxis.pdf" , bbox_inches='tight')
            else:
                fig.savefig( fig_output_dir + "/Patients/" +  patient +  ".pdf" , bbox_inches='tight')
            
            plt.close()
                


def prepare_data_rejection_analysis(fit_versions, res_dfs, min_DaysPostTransplant_forRejectionAnalysis, set_type):
    """
    preparing the matrices for ROC curve analysis
    """
    res_roc_dfs = res_dfs

    for fit_ver in fit_versions:
      
        res_df = res_roc_dfs[fit_ver].copy()
        
        if set_type == "lt":
            
            # multiplying by two single lung to compare cfDNA level
            res_df.ix[:,'DonorFraction_G1'] = res_df['DonorFraction_G1'] * (3 - res_df['Num_of_Lungs'])
            res_df.ix[:,'DonorFraction_G2'] = res_df['DonorFraction_G2'] * (3 - res_df['Num_of_Lungs'])
            res_df['DonorFraction_G2_CorrectedByError'] = res_df['DonorFraction_G2_CorrectedByError'] * (3 - res_df['Num_of_Lungs'])
            
            
            # filtering the relevant samples for the analysis
            res_df = res_df.ix[res_df.DaysPostTransplant >= min_DaysPostTransplant_forRejectionAnalysis,]
            res_df = res_df.ix[res_df.Biopsy_level != "No data",]
            
        elif set_type == "ht": 
            # filtering the relevant samples for the analysis
            res_df = res_df.ix[res_df.DaysPostTransplant >= min_DaysPostTransplant_forRejectionAnalysis,]
            res_df = res_df.ix[res_df.Biopsy_level != "No data",]
            
        elif set_type == "bm":
            pass
        else:
            raise ValueError("unknown set name:" + set_type)
        
    
        res_roc_dfs[fit_ver] = res_df
    
    return(res_roc_dfs)
    
    
def plot_rejection_AUC_comparison(fit_versions, res_roc_dfs, fig_output_dir, set_type, use_correctedG2 = False):
    """
    drawing ROC for ht and lt
    """

    for fit_ver in fit_versions:
        
        res_df = res_roc_dfs[fit_ver]
        
        if set_type == "lt":
            roc_pl_df = pd.DataFrame({"Method" : ["Two-genomes", "One genome"] * 4,
                              "Comparison" : "", "AUC" : np.nan, 
                              "AUC_down" : np.nan, "AUC_up" : np.nan, 
                              "pvalue" : np.nan})       
            
            # 1
            biopsy_level_all_labels = ["1. Normal","4. Moderate","5. Severe"] 
            biopsy_level_true_labels = ["4. Moderate","5. Severe"]
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Biopsy 0 vs. 3-4", ind = 0, use_correctedG2 = use_correctedG2)
    
            # 2
            biopsy_level_all_labels = ["1. Normal","3. Mild","4. Moderate","5. Severe"] 
            biopsy_level_true_labels = ["3. Mild","4. Moderate","5. Severe"]
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Biopsy 0 vs. 2-4", ind = 2, use_correctedG2 = use_correctedG2)
    
            # 3
            biopsy_level_all_labels = ["1. Normal","2. Minimal","3. Mild","4. Moderate","5. Severe"] 
            biopsy_level_true_labels = ["3. Mild","2. Minimal","4. Moderate","5. Severe"]
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Biopsy 0 vs. 1-4", ind = 4, use_correctedG2 = use_correctedG2)
    
            # 4 
            response = (res_df.AMR[res_df.AMR.isin([0.0,1.0])].isin([1.0])*1).tolist()
            pred_g1 = res_df.DonorFraction_G1[res_df.AMR.isin([0.0,1.0])].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.AMR.isin([0.0,1.0])].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.AMR.isin([0.0,1.0])].tolist()
    
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "AMR negative vs. positive", ind = 6, use_correctedG2 = use_correctedG2)
    
        elif set_type == "ht": 
            
            roc_pl_df = pd.DataFrame({"Method" : ["Two-genomes", "One genome"] * 4,
                              "Comparison" : "", "AUC" : np.nan, 
                              "AUC_down" : np.nan, "AUC_up" : np.nan, 
                              "pvalue" : np.nan})
    
            # 1
            biopsy_level_all_labels = ['1. Quiescence', '3. Moderate or severe rejection'] 
            biopsy_level_true_labels = ['3. Moderate or severe rejection']
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df,
                                           comparison = "Quiescence vs. Moderate or severe", ind = 0, use_correctedG2 = use_correctedG2)
    
            # 2
            biopsy_level_all_labels = ['2. Mild rejection', '3. Moderate or severe rejection'] 
            biopsy_level_true_labels = ['3. Moderate or severe rejection']
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Mild vs. Moderate or severe", ind = 2, use_correctedG2 = use_correctedG2)
    
            # 3
            biopsy_level_all_labels = ['1. Quiescence','2. Mild rejection'] 
            biopsy_level_true_labels = ['2. Mild rejection']
    
            response = (res_df.Biopsy_level[res_df.Biopsy_level.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Quiescence vs. Mild", ind = 4, use_correctedG2 = use_correctedG2)
    
            # 4
            res_df['Biopsy_level_roc'] = res_df['Biopsy_level']
            res_df.ix[res_df.Is_severe_rejection == True,'Biopsy_level_roc'] = '4. Severe rejection'
    
            biopsy_level_all_labels = ['1. Quiescence', '4. Severe rejection'] 
            biopsy_level_true_labels = ['4. Severe rejection']
    
            response = (res_df.Biopsy_level_roc[res_df.Biopsy_level_roc.isin(biopsy_level_all_labels)].isin(biopsy_level_true_labels)*1).tolist()            
            pred_g1 = res_df.DonorFraction_G1[res_df.Biopsy_level_roc.isin(biopsy_level_all_labels)].tolist()
            pred_g2 = res_df.DonorFraction_G2[res_df.Biopsy_level_roc.isin(biopsy_level_all_labels)].tolist()
            pred_g2cor = res_df.DonorFraction_G2_CorrectedByError[res_df.Biopsy_level_roc.isin(biopsy_level_all_labels)].tolist()
    
            roc_pl_df = cal_roc_and_set_df(response, pred_g1, pred_g2,  pred_g2cor, roc_pl_df, 
                                           comparison = "Quiescence vs. Severe", ind = 6, use_correctedG2 = use_correctedG2)
    
        elif set_type == "bm":
            pass
        else:
            raise ValueError("unknown set name:" + set_type)
    
    
        if set_type in ['lt','ht']:
            if (use_correctedG2):
                r_plot_auc(roc_pl_df ,robjects.StrVector([fig_output_dir + "/" + fit_ver + "_AUC_bar_corrected" + ".pdf"]))
                roc_pl_df.to_csv(fig_output_dir + "/" + fit_ver + "_AUC_corrected" + ".tsv", sep='\t', index = False)
                
            else:
                r_plot_auc(roc_pl_df ,robjects.StrVector([fig_output_dir + "/" + fit_ver + "_AUC_bar_noncorrected" + ".pdf"]))
                roc_pl_df.to_csv(fig_output_dir + "/" + fit_ver + "_AUC_noncorrected" + ".tsv", sep='\t', index = False)
    


###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Comparing inference of donor-derived cfDNA fraction using recipient and donor genotypes and only recipient")
    
    parser.add_argument("output_dir", type=str, help= "The output directory of the analysis, this directory should contain the fit directory")

    parser.add_argument("set_type", type=str, help="set_type: lt,ht or bm")
    
    parser.add_argument("fig_output_dir", type=str, help= "the name of the directory inside the analysis directory Output/Analysis/<dir>")

    parser.add_argument("-f", "--fit_versions", dest='fit_versions', default=['V21', 'V22', 'V23', 'V24'], nargs='*',
                      help=('fit verions (for example: V21), can be multiple strings'))
    
    parser.add_argument("-d", "--model_description_file", dest='model_description_file', default=None, type=str,
                      help=('filename with model names'))
    
    parser.add_argument("-p", "--draw_patients", dest='draw_patients',  action='store_true',
                        help=('draw patients plots'))
    
    
    parser.add_argument("-b", "--fix_by_BC_flag", dest='fix_by_BC_flag',  action='store_true',
                        help=('for each patient correct predictions of dd-cfDNA in sample before the transplant'))
    
    
    
    
    
    
    
    args = parser.parse_args()
    
    output_dir = args.output_dir
    set_type = args.set_type
    fit_versions = args.fit_versions
    fig_output_dir = output_dir + '/Analysis/' + args.fig_output_dir

    if not os.path.exists(fig_output_dir):
        os.makedirs(fig_output_dir)

    if not os.path.exists(fig_output_dir + '/Patients'):
        os.makedirs(fig_output_dir + '/Patients')


    ###############################################################################################################
    # parameters (TODO - make options of the scripts? these are fixed to be similar to previous 2G papers)
    ###############################################################################################################
    
    
    
    if set_type == "lt":
        min_2G_seq_error = 0.5
        min_DaysPostTransplant = -10
        min_reads_that_map_to_snps = 1e5
        min_DaysPostTransplant_forRejectionAnalysis = 60
    elif set_type == "ht":
        min_2G_seq_error = 0.15
        min_DaysPostTransplant = 15
        min_reads_that_map_to_snps = 1e5
        min_DaysPostTransplant_forRejectionAnalysis = 15
    elif set_type == "bm":
        min_2G_seq_error = 1
        min_DaysPostTransplant = -100
        min_reads_that_map_to_snps = 1e5
        min_DaysPostTransplant_forRejectionAnalysis = -100
    else:
        raise ValueError("unknown set name:" + set_type)


    print("Running comparison version 0.1")
    
    models_desc_dict = {}
    for model in args.fit_versions:
        models_desc_dict[model] = model
    
    if args.model_description_file is not None:
        print("Loading model description")
        model_desc_df = pd.read_table(args.model_description_file, sep='\t', na_values = "")
        for r,row in model_desc_df.iterrows():
            print(row['Model'])
            models_desc_dict[row['Model']] = row['Model_description']
            
        
    print(models_desc_dict)  
    
    print("Loading the data")
    res_dfs = load_data(output_dir, fit_versions, corrupted_sampels, min_2G_seq_error, min_DaysPostTransplant, min_reads_that_map_to_snps)
    
    print("Parsing clinical terms")
    res_dfs = parse_clinical_terms(fit_versions, res_dfs, set_type)
    
    print("Fix 2G by sequencing error")
    res_dfs = fix_2g_by_seq_error(fit_versions, res_dfs, set_type, args.fix_by_BC_flag)
    
    print("Plot 2G vs. 1G")
    plot_2G_vs_1G(fit_versions, res_dfs, fig_output_dir, set_type)
    
    print("Preparing data for rejection analysis")
    res_roc_dfs = prepare_data_rejection_analysis(fit_versions, res_dfs, min_DaysPostTransplant_forRejectionAnalysis, set_type)
    
    print("Ploting rejection AUC comparison (2G not corrected)")
    plot_rejection_AUC_comparison(fit_versions, res_roc_dfs, fig_output_dir, set_type, use_correctedG2 = False)
    plot_rejection_AUC_comparison(fit_versions, res_roc_dfs, fig_output_dir, set_type, use_correctedG2 = True)
    
    print("Compare 1G algorithms (i.e. parameters)")
    plot_2G_vs_1G_algorithms(fit_versions, res_dfs, fig_output_dir, set_type, models_desc_dict)
    
    if (args.draw_patients):
        print("Patients plots no axes labels")
        plot_patients(fit_versions, res_dfs, fig_output_dir, set_type, models_desc_dict,use_axis_labels = False)
        print("Patients plots with axes labels")
        plot_patients(fit_versions, res_dfs, fig_output_dir, set_type, models_desc_dict,use_axis_labels = True)
        
    
    
    print('Done')

    
if __name__ == '__main__':
  main()



