#!/usr/bin/env python3

"""
Drawing simulations results
Eilon Sharon 2017
"""

#~/postdoc/cfDNA/bm_sim20 python ~/bin/cfDNA1G/simulations/cfDNA_draw_sim_results.py ./bm.samples.tsv Output Vmulti -f V37 -a ../bm_sim15/Output ../bm_sim20HaoWithoutU/Output ../bm_sim21/Output ../bm_sim22/Output

import os
import argparse
import gzip

import itertools

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy import interpolate


import seaborn as sns
sns.set(rc={"figure.figsize": (6, 6)})

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

###############################################################################################################
# constants
###############################################################################################################


#cbPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# The palette with grey:
#cbgPalette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# The palette with black:
#cbbPalette = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# The palette with grey:
cbgPalette = ["#999999", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9"]
# The palette with black:
cbbPalette = ["#000000", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9"]

cbPalette = [ "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#56B4E9"]

###############################################################################################################
# helper functions
###############################################################################################################

def load_simulation_results(org_samples_table_filename, output_dir, fit_versions, additional_output_dirs):


    print("Loading original samples table")
    samples_df_filename = os.path.expanduser(org_samples_table_filename)
    org_samples_df = pd.read_table(samples_df_filename, sep='\t', na_values = "") 
    org_samples_df['DonorCorrected'] = org_samples_df['Donor'] / 100.0

    res_dfs = {}
    
    for fit_ver in fit_versions:
      
      
      
        res_df = load_simulation_dataframe(output_dir + '/Fit/fit_sim_' + fit_ver + '.res.tsv', fit_ver,
                                           org_samples_df)
        
        for add_output in additional_output_dirs:
          cur_input_filename = add_output + '/Fit/fit_sim_' + fit_ver + '.res.tsv'
          if (os.path.exists(cur_input_filename)):
            add_res_df = load_simulation_dataframe(cur_input_filename, fit_ver, 
                                                 org_samples_df)
            res_df = pd.concat([res_df, add_res_df])
          else:
            print("File does not exists: %s" % (cur_input_filename))

        #TODO - remove. this fixes a bug in the input tables
        res_df.ix[res_df.Sim_NumSNPs > 2000000,'Sim_NumSNPs'] = 2000000

        res_dfs[fit_ver] = res_df
        
        
        

        
        
    return(samples_df_filename, res_dfs)
  

def load_simulation_dataframe(input_filename, fit_ver, org_samples_df):
  
  res_df = pd.read_table(input_filename, sep='\t', na_values = "")   
  print("%s number of samples: %d (file: %s)" %(fit_ver, res_df.shape[0], input_filename))
  
  samples_2g_df  = org_samples_df[ ['SampleUniqueStr','DonorCorrected'] ].copy() 
  samples_2g_df.ix[samples_2g_df.SampleUniqueStr.str.contains('_BC'),'DonorCorrected'] = 0.0
  samples_2g_df = samples_2g_df.rename(index=str, columns={"SampleUniqueStr" : "s1_SampleUniqueStr",  "DonorCorrected": "Sample1_dd_cfDNA_G2"})
  res_df = res_df.merge(samples_2g_df, how='left', left_on='Sample1', right_on='s1_SampleUniqueStr')
  samples_2g_df = samples_2g_df.rename(index=str, columns={"s1_SampleUniqueStr" : "s2_SampleUniqueStr","Sample1_dd_cfDNA_G2": "Sample2_dd_cfDNA_G2"})
  res_df = res_df.merge(samples_2g_df, how='left', left_on='Sample2', right_on='s2_SampleUniqueStr')
  
  res_df['expected_dd_cfDNA'] = ( res_df['Sample1_dd_cfDNA_G2'] * (1.0 - res_df['Sample2_fraction']) +
                                  res_df['Sample2_dd_cfDNA_G2'] * (      res_df['Sample2_fraction']) )

  # if closely related then assuming that two pure (_BC) samples
  res_df.ix[ (res_df['Sample1'].str.contains('_BC') & (res_df['Sample2'].str.contains('_BC')) )  ,'expected_dd_cfDNA'] = res_df['Sample2_fraction']

  print('Results exists for %d samples' % (np.sum(~res_df.OneGenome_P_donor.isnull())))
  return(res_df)


def plot_2_fit_ver_comparison(res_df1, res_df2, fit_ver1, fit_ver2, output_plot_filename):
  """
  comparing to fit versions
  """
  
  res_df1 = res_df1.copy()
  res_df2 = res_df2.copy()
  
  res_df2 = res_df2.ix[:, ['SampleUniqueStr', 'OneGenome_P_donor'] ]
  
  
  res_df1 = res_df1.rename(index=str, columns={"OneGenome_P_donor": "OneGenome_P_donor1"})
  res_df2 = res_df2.rename(index=str, columns={"OneGenome_P_donor": "OneGenome_P_donor2"})

  res_df = res_df1.merge(res_df2, on = 'SampleUniqueStr', how='inner')
  
  
  cur_x_values = res_df.OneGenome_P_donor1.values    * 100
  cur_y_values = res_df.OneGenome_P_donor2.values * 100
                 
  
  max_val = max(np.max(cur_x_values), np.max(cur_y_values))
  
  with sns.axes_style("white"):
    
    sns.set_style("ticks")
    fig, ax = plt.subplots(1,1,figsize=(6,6))
  
    ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None, linewidth = 1.0)
    
    colors_vec = ["#001f3f","#39CCCC","#3D9970","#2ECC40","#FF851B","#0074D9","#FF4136","#85144b","#F012BE","#B10DC9","#AAAAAA","#01FF70","#FFDC00","#111111","#DDDDDD"]
    
    colors = itertools.cycle(sns.color_palette(colors_vec))
    cur_color = next(colors)
    
    ax.plot(cur_x_values, cur_y_values, 
            label=None, marker='o', linestyle='', color=cur_color,
            markersize=3,
            fillstyle='full',
            alpha=1)
    
    plt.xscale('log')#, nonposy='clip')
    plt.yscale('log')#, nonposy='clip')
    
    
    
    plt.xlabel(fit_ver1 + ": " + "Predicted percent dd-cfDNA (log scale)",fontsize=12)
    plt.ylabel(fit_ver2 + ": " + "Predicted percent dd-cfDNA (log scale)",fontsize=12)
    
    spearman_rho, _ = spearmanr(cur_x_values, cur_y_values)
    pearson_r, _ = pearsonr(cur_x_values, cur_y_values)
    
    
    ax.text(0.9, 0.02, r"Spearman's $\rho$  = %.4f" % (spearman_rho),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=10)

    fig.savefig(output_plot_filename , bbox_inches='tight')
    plt.close()
    
  #print(res_df.ix[res_df.OneGenome_P_donor2 < 1e-3 ,["SampleUniqueStr", "OneGenome_P_donor1", "OneGenome_P_donor2", "expected_dd_cfDNA"]])
    
  
  res_df.to_csv(output_plot_filename.replace(".pdf",".tsv"), sep='\t', index = False)
    
        
  

def plot_sim_results(res_df, output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, draw_median_line = True, use_markersAndFill_for_snpXread= False, use_markersAndFill_for_relatedness= False):
    """
    type_of_y_axis = prediction | mae | cdf
    """
    
    #colors = itertools.cycle(sns.color_palette("cubehelix", len(reads_per_snp))[::-1])
    #marker_fill_styles = itertools.cycle(['left','right','bottom','top','none'])
    markers_vec = ["o","s","d","X","8","*","h",">","v","<","^","h","P", "H", "D" ] #('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
    #markers = itertools.cycle(markers_vec)
    
    #"#7FDBFF",
    colors_vec = ["#001f3f","#39CCCC","#3D9970","#2ECC40","#FF851B","#0074D9","#FF4136","#85144b","#F012BE","#B10DC9","#AAAAAA","#01FF70","#FFDC00","#111111","#DDDDDD"]
    colors_vec_perm = ["#001f3f","#85144b","#2ECC40","#B10DC9","#FF851B","#3D9970","#0074D9","#39CCCC","#FF4136","#F012BE","#AAAAAA","#01FF70","#FFDC00","#111111","#DDDDDD"]
    linestyle_vec = ['-', '--', '-.', ':']
    
    res_df = res_df.copy()
    res_df.reset_index(inplace=True)
    
    reads_per_snp = res_df.Sim_reads_per_SNP.unique()
    reads_per_snp.sort()
    reads_per_snp = reads_per_snp[::-1]
    
    total_snp_num = res_df.Sim_NumSNPs.astype(int).unique() # OneGenome_number_of_SNPs
    total_snp_num.sort()
    total_snp_num = total_snp_num[::-1]
    
    res_df['relatedness'] = "Separated by: " + res_df.ibd_m1.astype(str) + " and " + res_df.ibd_m1.astype(str) + "meioses events"
    res_df.ix[ (res_df.ibd_m1 == 100) & (res_df.ibd_m2 == 100)   ,'relatedness'] = "1. unrelated"
    res_df.ix[ (res_df.ibd_m1 == 2) & (res_df.ibd_m2 == 2)   ,'relatedness'] = "2. siblings"
    res_df.ix[ (res_df.ibd_m1 == 2) & (res_df.ibd_m2 == 100)   ,'relatedness'] = "3. half-siblings"
    res_df.ix[ (res_df.ibd_m1 == 4) & (res_df.ibd_m2 == 100)   ,'relatedness'] = "4. first cousins"
    res_df.ix[ (res_df.ibd_m1 == 6) & (res_df.ibd_m2 == 100)   ,'relatedness'] = "5. second cousins"
    
    
    res_df['plotlabel'] = "non"
    
    relatedness_types = res_df.relatedness.unique()
    relatedness_types.sort()
    
    
    
    max_val = max(res_df.expected_dd_cfDNA.max(), res_df.OneGenome_P_donor.max()) * 100
    spearman_rho, _ = spearmanr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    pearson_r, _ = pearsonr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    
    print("In plot %s:" % (output_plot_filename))
    print("MAE (mean): %f" % ( np.mean(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)*100) ))
    print("MAE (median): %f" % ( np.median(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)*100) ))
    
    print("spearman rho: %f" % (spearman_rho)) 
    print("pearson R: %f" % (pearson_r)) 
    
    

    
    with sns.axes_style("white"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(8,8))
        
        if (type_of_y_axis == "prediction"):
          ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None, linewidth = 1.0)
        elif (type_of_y_axis == "mae"):
          ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0,max_val*0.3], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0,max_val*0.1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0,max_val*0.05], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0,max_val*0.01], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0.01,0.01], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0,max_val],[0.5,0.5], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0,max_val],[0.1,0.1], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0,max_val],[1,1], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0,max_val],[5,5], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0,max_val],[10,10], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0.1,0.1],[0,100], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([0.5,0.5],[0,100], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([1,1],[0,100], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([5,5],[0,100], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          ax.plot([10,10],[0,100], linestyle=':', color = "gray", label=None, linewidth = 0.5)
          
        elif (type_of_y_axis == "cdf"):
          ax.plot([0,max_val],[0.9,0.9], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0.8,0.8], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0.5,0.5], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.001,0.001],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.01,0.01],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.1,0.1],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.2,0.2],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([1,1],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
        elif (type_of_y_axis == "relcdf"):
          ax.plot([0,max_val],[0.9,0.9], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0.8,0.8], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0,max_val],[0.5,0.5], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.001,0.001],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.01,0.01],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.1,0.1],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([0.5,0.5],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          ax.plot([1,1],[0,1], linestyle=':', color = "gray", label=None, linewidth = 1.0)
          
          
        
        colors_for_all = itertools.cycle(sns.color_palette(colors_vec)) # "hls",10 #sns.color_palette("cubehelix", len(reads_per_snp))[::-1]
            
        colors = itertools.cycle(colors_vec_perm) #(cbPalette) #sns.color_palette("cubehelix", len(reads_per_snp))[::-1]
        marker_fill_styles = itertools.cycle(['left','right','bottom','top'])
        markers = itertools.cycle(markers_vec)
        line_styles = itertools.cycle(linestyle_vec)
        
        for cur_relatedness in  relatedness_types:
          
          if (not use_markersAndFill_for_snpXread):
            
            if (not use_markersAndFill_for_relatedness):
              if (only_unrelated):
                cur_marker = 'o'
              else:
                cur_color = next(colors)
                
              markers = itertools.cycle(markers_vec)
              line_styles = itertools.cycle(linestyle_vec)
            
            else:
              cur_marker_fillstyle = 'full'
              cur_color = next(colors)
              cur_marker = next(markers)
              
              
            
          for design_num_snps in total_snp_num:
            
            if ((not use_markersAndFill_for_snpXread) and (not use_markersAndFill_for_snpXread)):
              
              if (only_unrelated):
                cur_marker_fillstyle = next(marker_fill_styles)
              else:
                cur_marker = next(markers)
                
              cur_line_style = next(line_styles)
              
              if (only_unrelated):
                colors = itertools.cycle(colors_vec_perm) #(cbPalette) #sns.color_palette("cubehelix", len(reads_per_snp))[::-1]
              else:
                marker_fill_styles = itertools.cycle(['left','right','bottom','top'])
              
            for r in (reads_per_snp):
                
                if ((not use_markersAndFill_for_snpXread) and (not use_markersAndFill_for_snpXread)):
                  
                  if (only_unrelated):
                    cur_color = next(colors)
                  else:
                    cur_marker_fillstyle = next(marker_fill_styles)
                elif (use_markersAndFill_for_snpXread):
                  cur_color = next(colors)
                  cur_marker = next(markers)
                  cur_marker_fillstyle = 'full'
                  
                  
                
                cur_markeredgecolor = None
                cur_markeredgewidth = None
                
                #if r == 2.0:
                #    cur_markeredgewidth = 1
                #    cur_markeredgecolor = 'black'
                #    cur_marker_fill_style = 'none'
                    
                
                cur_inds_tf = ((res_df.Sim_reads_per_SNP == r) &
                               (res_df.Sim_NumSNPs == design_num_snps) &
                               (res_df.relatedness == cur_relatedness))
                
                #print(np.nonzero(cur_inds_tf))
                
                if (np.sum(cur_inds_tf) > 0):
                
                  cur_x_values = res_df.expected_dd_cfDNA[cur_inds_tf].values * 100
                  cur_y_values = res_df.OneGenome_P_donor[cur_inds_tf].values * 100
                  
                  if (type_of_y_axis == "mae"):
                    cur_y_values = np.abs(cur_y_values - cur_x_values)
                  elif (type_of_y_axis == "cdf"):
                    abs_error = np.sort(np.absolute(cur_y_values-cur_x_values))
                    abs_error_cdf = np.array(range(abs_error.shape[0])) * (1.0 / abs_error.shape[0])
                    cur_x_values = abs_error
                    cur_y_values = abs_error_cdf
                  elif (type_of_y_axis == "relcdf"):
                    abs_error = np.sort(np.absolute( (cur_y_values-cur_x_values)/np.maximum(cur_x_values,1e-6) ))
                    abs_error_cdf = np.array(range(abs_error.shape[0])) * (1.0 / abs_error.shape[0])
                    cur_x_values = abs_error
                    cur_y_values = abs_error_cdf
                    
                  cur_mean_num_snps = int(res_df.OneGenome_number_of_SNPs[cur_inds_tf].mean())
                  cur_mean_cov =  np.mean(res_df.OneGenome_number_of_reads_mapped_to_SNPs[cur_inds_tf] / res_df.OneGenome_number_of_SNPs[cur_inds_tf])
                  

                  cur_label = ("%d SNPs, %.2f reads per SNP" % (cur_mean_num_snps, cur_mean_cov))
                  if (not only_unrelated):
                    cur_label = cur_label + (", relatedness: %s" % (cur_relatedness))
                  
                  #print(np.flatnonzero(cur_inds_tf))
                  res_df.ix[np.flatnonzero(cur_inds_tf),'plotlabel'] = cur_label
                  
                  if (type_of_y_axis == "prediction"):
                    cur_spearman_rho, _ = spearmanr(cur_x_values,cur_y_values)
                    cur_label = cur_label + (r" (Spearman's $\rho$  = %.3f)" % (cur_spearman_rho) )


                  
                  
                  #print("ploting %s with %s,%s,%s" % (cur_label, cur_marker, cur_marker_fillstyle, str(cur_color)))
                  
                  #print(cur_label)
                  #print('Y' * 100)
                  #print(cur_y_values)
                  #print('X' * 100)
                  #print(cur_x_values)
                  
                  if (type_of_y_axis in ["cdf", "relcdf"]):
                    ax.plot(cur_x_values, cur_y_values, 
                          label=cur_label,  linestyle=cur_line_style, color=next(colors_for_all))
                  else:
                    ax.plot(cur_x_values, np.maximum(6*1e-5,cur_y_values), 
                          label=None, marker=cur_marker, linestyle='', color=cur_color,
                          markersize=3,
                          fillstyle=cur_marker_fillstyle,
                          markeredgecolor = cur_markeredgecolor,
                          markeredgewidth = cur_markeredgewidth,
                          alpha=0.6)
                    
                    # drawing mean line
                    cur_x_sort_inds = np.argsort(cur_x_values)
                    #print('Z' * 100)
                    
                    cur_mean_df = pd.DataFrame({"x" : cur_x_values[cur_x_sort_inds], "y" : cur_y_values[cur_x_sort_inds]})
                    #print(cur_mean_df)
                    cur_mean_df = cur_mean_df.groupby("x",  as_index = False).median()
                    
                    #print('D' * 100)
                    #print(cur_mean_df.y.values)
                    if (draw_median_line):
                      cur_line_style = '-'
                    else:
                      cur_line_style = ''
                      
                    #ax.plot(cur_mean_df.x.values, cur_mean_df.y.values, 
                    #      label=cur_label, marker=cur_marker, linestyle=cur_line_style, color=cur_color, linewidth = 0.5,
                    #      markersize=5,
                    #      fillstyle=cur_marker_fillstyle,
                    #      markeredgecolor = cur_markeredgecolor,
                    #      markeredgewidth = cur_markeredgewidth,
                    #      alpha=1.0)
                    
                    #cur_mean_df = cur_mean_df.ix[(cur_mean_df.x>0) & (cur_mean_df.y>0),:]
                    inter_x = cur_mean_df.x.unique()
                    
                    #print(cur_x_values[cur_x_sort_inds])
                    #print(cur_y_values[cur_x_sort_inds])
                    #print(len(inter_x))# - np.sqrt(2*len(inter_x)))
                    #ius = interpolate.InterpolatedUnivariateSpline(cur_x_values[cur_x_sort_inds], cur_y_values[cur_x_sort_inds]) #s=len(inter_x) - 2*np.sqrt(len(inter_x))  , s=0
                    cur_s = np.sum((cur_mean_df.x.values > 0) & (cur_mean_df.y.values > 0))
                    cur_s = cur_s / 5
                    #print(cur_s)
                    
                    cur_spline = interpolate.UnivariateSpline(x=np.log10(cur_mean_df.x.values[(cur_mean_df.x.values > 0) & (cur_mean_df.y.values > 1e-4)]), y=np.log10(cur_mean_df.y.values[(cur_mean_df.x.values > 0) & (cur_mean_df.y.values > 1e-4) ]), s=cur_s) # , s=0
                    #print(cur_spline)
                    
                    #print(np.log10(inter_x))
                    inter_y = 10**cur_spline(np.log10(inter_x))
                    #print(inter_y)
                    #print(len(cur_x_values[cur_x_sort_inds]))
                    #print(len(inter_x))
                    ax.plot(cur_mean_df.x.values, np.maximum(6*1e-5,cur_mean_df.y.values), 
                          label=cur_label, marker=cur_marker, linestyle='', color=cur_color, linewidth = 0.5,
                          markersize=5,
                          fillstyle=cur_marker_fillstyle,
                          markeredgecolor = cur_markeredgecolor,
                          markeredgewidth = cur_markeredgewidth,
                          alpha=1.0)
                    
                    if (type_of_y_axis == "mae"):
                      ax.plot(inter_x, inter_y, 
                            label=None, marker='', linestyle=cur_line_style, color=cur_color, linewidth = 1.0,
                            markersize=5,
                            fillstyle=cur_marker_fillstyle,
                            markeredgecolor = cur_markeredgecolor,
                            markeredgewidth = cur_markeredgewidth,
                            alpha=1.0)
                      
                      print(cur_label)
                      print(np.mean(cur_y_values))
                      print(np.median(cur_y_values))
                    
                    #fit = np.polyfit(np.log10(cur_x_values[cur_x_sort_inds]), np.log10(cur_y_values[cur_x_sort_inds]), deg=1)
                    
                    #ax.plot(cur_mean_df.x.values, 10**(fit[0] * np.log10(cur_mean_df.x.values) + fit[1]), 
                    #      label=cur_label, marker='', linestyle=cur_line_style, color=cur_color, linewidth = 1,
                    #      markersize=5,
                    #      fillstyle=cur_marker_fillstyle,
                    #      markeredgecolor = cur_markeredgecolor,
                    #      markeredgewidth = cur_markeredgewidth,
                    #      alpha=1.0)
                    
                  
                  
                    
                
                  
                  if (type_of_y_axis == "mae"):
                    cur_inds = np.argsort(cur_x_values)
                    cur_x_values = cur_x_values[cur_inds]
                    cur_y_values = cur_y_values[cur_inds]
                    
                    cur_smoothed_y = np.convolve(cur_y_values, np.ones(5)/5, 'same')
                    
                    #ax.plot(cur_x_values, cur_smoothed_y, 
                    #      label=None, marker=cur_marker, linestyle='-', color=cur_color,
                    #      markersize=6,
                    #      fillstyle=cur_marker_fillstyle,
                    #      markeredgecolor = cur_markeredgecolor,
                    #      markeredgewidth = cur_markeredgewidth,
                    #      alpha=0.7)
        
        if (type_of_y_axis in ["cdf", "relcdf"]):
          plt.xlim(1e-6, plt.xlim()[1])
        else:
          plt.ylim(5*1e-5, 1e2)
        
        if (log_scale):
          
          #if (type_of_y_axis != "relcdf"):
          plt.xscale('log', nonposy='clip')
          if (type_of_y_axis not in ["cdf","relcdf"]):
            plt.yscale('log', nonposy='clip')
          
        plt.xlabel("Simulated percent dd-cfDNA (log scale)",fontsize=12)
        if (type_of_y_axis == "prediction"):
          plt.ylabel("Predicted percent dd-cfDNA (log scale)",fontsize=12)
          plt.legend(title="",loc='lower right',fontsize=8) 
        elif (type_of_y_axis == "mae"):  
          plt.ylabel("Absolute error of predicted percent dd-cfDNA (log scale)",fontsize=12)
          plt.legend(title="",loc='upper left',fontsize=8)
        elif (type_of_y_axis == "cdf"):  
          plt.xlabel("Absolute error of predicted percent dd-cfDNA (log scale)",fontsize=12)
          plt.ylabel("Samples with smaller or equal absolute error (fraction)",fontsize=12)
          plt.legend(title="",loc='upper left',fontsize=7)
        elif (type_of_y_axis == "relcdf"):  
          plt.xlabel("Relative absolute error of predicted precent dd-cfDNA (log scale)",fontsize=12)
          plt.ylabel("Samples with smaller or equal relative absolute error (fraction)",fontsize=12)
          plt.legend(title="",loc='upper left',fontsize=7)
        
        sns.despine()
        
        #cur_ticker = matplotlib.ticker.FuncFormatter(lambda x,pos: '%f' % (x)  )
        
        #ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
         
        #ax.text(0.1, 0.98, r"Pearson $R^2$ = %.2f\\Spearman's $\rho$  = %.2f" % (pearson_r**2, spearman_rho),
        #        verticalalignment='bottom', horizontalalignment='right',
        #        transform=ax.transAxes,
        #        color='black', fontsize=10)
        
        fig.savefig( output_plot_filename , bbox_inches='tight')
        plt.close()
        
        
        #fig, ax = plt.subplots(1,1,figsize=(8,8))
        
        #res_df['expectedddcfDNA100'] = res_df.expected_dd_cfDNA*100.0
        #res_df['OneGenomePdonor100'] = res_df.OneGenome_P_donor*100.0
        #res_df['unit'] = np.linspace(0,res_df.shape[0]-1,res_df.shape[0])
        
        #print(res_df['plotlabel'] )
        
        #sns.tsplot(res_df, time="expectedddcfDNA100",  unit="unit", condition="plotlabel", value="OneGenomePdonor100", err_style='ci_band', ci=95, legend=True, ax=ax) #interpolate=True, color=None, estimator=<function mean>, n_boot=5000, err_palette=None, err_kws=None, 
        #if (log_scale):
        #  plt.xscale('log', nonposy='clip')
        #  plt.yscale('log', nonposy='clip')
        #fig.savefig( output_plot_filename.replace(".pdf","_smooth.pdf") , bbox_inches='tight')
        #plt.close()
        
        

 
  
def plot_sim_results_dotplot(res_df, output_plot_filename, log_scale = True):
    
    reads_per_snp = res_df.Sim_reads_per_SNP.unique()
    reads_per_snp.sort()
    max_val = max(res_df.expected_dd_cfDNA.max()*100, res_df.OneGenome_P_donor.max()*100)
    spearman_rho, _ = spearmanr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    pearson_r, _ = pearsonr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    
    print("In plot %s:" % (output_plot_filename))
    print("MAE (mean): %f" % ( np.mean(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)) ))
    print("MAE (median): %f" % ( np.median(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)) ))
    
    print("spearman rho: %f" % (spearman_rho)) 
    print("pearson R: %f" % (pearson_r)) 
    
    colors = itertools.cycle(sns.color_palette("cubehelix", len(reads_per_snp))[::-1])
    marker_fill_styles = itertools.cycle(['left','right','bottom','top'])
    markers_vec = ["o","s","d","+","x","8","h",">","v","<","^","h" ] 
    markers = itertools.cycle(markers_vec)
    

    
    with sns.axes_style("white"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(6,6)) 
        ax.plot([0,max_val],[0,max_val], linestyle=':', color = "gray", label=None)
        handles = [None] * len(reads_per_snp)
        for r in (reads_per_snp):
            cur_color = next(colors)
            cur_marker_fill_style = next(marker_fill_styles)
            cur_markeredgecolor = None
            cur_markeredgewidth = None
            if r == 2.0:
                cur_markeredgewidth = 1
                cur_markeredgecolor = 'black'
                cur_marker_fill_style = 'none'
            ax.plot(res_df.expected_dd_cfDNA[res_df.Sim_reads_per_SNP == r]*100, res_df.OneGenome_P_donor[res_df.Sim_reads_per_SNP == r]*100, 
                    label=str(r), marker='o', linestyle='', color=cur_color, markersize=6, fillstyle=cur_marker_fill_style, markeredgecolor = cur_markeredgecolor, markeredgewidth = cur_markeredgewidth)
        if (log_scale):
          plt.xscale('log', nonposy='clip')
          plt.yscale('log', nonposy='clip')
        plt.xlabel("Simulated dd-cfDNA (log scale)",fontsize=12)
        plt.ylabel("Predicted dd-cfDNA (log scale)",fontsize=12)
        plt.legend(title="Avg. coverage",loc='upper left') 
        sns.despine()
        
        #cur_ticker = matplotlib.ticker.FuncFormatter(lambda x,pos: '%f' % (x)  )
        
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
         
        ax.text(0.9, 0.02, r"Pearson $R^2$ = %.2f\\Spearman's $\rho$  = %.2f" % (pearson_r**2, spearman_rho),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=10)
        
        fig.savefig( output_plot_filename , bbox_inches='tight')
        plt.close()
   
def plot_sim_results_mae(res_df, output_plot_filename):
    
    reads_per_snp = res_df.Sim_reads_per_SNP.unique()
    reads_per_snp.sort()
    max_val = max(res_df.expected_dd_cfDNA.max()*100, res_df.OneGenome_P_donor.max()*100)
    spearman_rho, _ = spearmanr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    pearson_r, _ = pearsonr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    
    print("In plot %s:" % (output_plot_filename))
    print("MAE (mean): %f" % ( np.mean(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)) ))
    print("MAE (median): %f" % ( np.median(np.absolute(res_df.expected_dd_cfDNA-res_df.OneGenome_P_donor)) ))
    
    print("spearman rho: %f" % (spearman_rho)) 
    print("pearson R: %f" % (pearson_r)) 
    
    colors = itertools.cycle(sns.color_palette("cubehelix", len(reads_per_snp))[::-1])
    marker_fill_styles = itertools.cycle(['left','right','bottom','top'])

    
    with sns.axes_style("white"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(6,6)) 
        ax.plot([0,max_val],[0.2,0.2], linestyle=':', color = "gray", label=None)
        ax.plot([0,max_val],[1.0,1.0], linestyle=':', color = "gray", label=None)
        ax.plot([0,max_val],[4.0,4.0], linestyle=':', color = "gray", label=None)
        handles = [None] * len(reads_per_snp)
        for r in (reads_per_snp):
            cur_color = next(colors)
            cur_marker_fill_style = next(marker_fill_styles)
            cur_markeredgecolor = None
            cur_markeredgewidth = None
            if r == 2.0:
                cur_markeredgewidth = 1
                cur_markeredgecolor = 'black'
                cur_marker_fill_style = 'none'
            ax.plot(res_df.expected_dd_cfDNA[res_df.Sim_reads_per_SNP == r]*100, 
                    np.abs((res_df.expected_dd_cfDNA[res_df.Sim_reads_per_SNP == r] - res_df.OneGenome_P_donor[res_df.Sim_reads_per_SNP == r])*100), 
                    label=str(r), marker='o', linestyle='', color=cur_color, markersize=6, fillstyle=cur_marker_fill_style, markeredgecolor = cur_markeredgecolor, markeredgewidth = cur_markeredgewidth)
        plt.xscale('log', nonposy='clip')
        plt.yscale('log', nonposy='clip')
        plt.xlabel("Simulated dd-cfDNA (log scale)",fontsize=12)
        plt.ylabel("Absolute error od dd-cfDNA prediction (log scale)",fontsize=12)
        plt.legend(title="Avg. coverage",loc='upper left') 
        sns.despine()
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
         
        ax.text(0.9, 0.02, r"Pearson $R^2$ = %.2f\\Spearman's $\rho$  = %.2f" % (pearson_r**2, spearman_rho),
                verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=10)
        
        fig.savefig( output_plot_filename , bbox_inches='tight')
        plt.close()
        

def plot_sim_results_cdf(res_df, output_plot_filename):
    
    reads_per_snp = res_df.Sim_reads_per_SNP.unique()
    reads_per_snp.sort()
    
    spearman_rho, _ = spearmanr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    pearson_r, _ = pearsonr(res_df.expected_dd_cfDNA, res_df.OneGenome_P_donor)
    

    
    colors = itertools.cycle(sns.color_palette("cubehelix", len(reads_per_snp))[::-1])
    marker_fill_styles = itertools.cycle(['left','right','bottom','top'])

    print("In plot %s:" % (output_plot_filename))
    
    with sns.axes_style("white"):
        
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(6,6)) 
   
        handles = [None] * len(reads_per_snp)
        for r in (reads_per_snp):
            cur_color = next(colors)
            cur_marker_fill_style = next(marker_fill_styles)
            cur_markeredgecolor = None
            cur_markeredgewidth = None
            if r == 2.0:
                cur_markeredgewidth = 1
                cur_markeredgecolor = 'black'
                cur_marker_fill_style = 'none'
                
            abs_error = np.sort(np.absolute(res_df.expected_dd_cfDNA[res_df.Sim_reads_per_SNP == r] - res_df.OneGenome_P_donor[res_df.Sim_reads_per_SNP == r])*100)
            abs_error_cdf = np.array(range(abs_error.shape[0])) * (1.0 / abs_error.shape[0])
            ax.plot(abs_error, abs_error_cdf,
                    label=str(r), marker='', linestyle='-', color=cur_color, markersize=6, fillstyle=cur_marker_fill_style, markeredgecolor = cur_markeredgecolor, markeredgewidth = cur_markeredgewidth)
        #plt.xscale('log', nonposy='clip')
        #plt.yscale('log', nonposy='clip')
        plt.xlabel("Absolute difference between the predicted and expected precent of dd-cfDNA(\%)",fontsize=12)
        plt.ylabel("Samples with smaller or equal absolute error",fontsize=12)
        plt.legend(title="Avg. coverage",loc='lower right') 
        sns.despine()
        plt.xlim((0,2.0))
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        fig.savefig( output_plot_filename , bbox_inches='tight')
        plt.close()
        
        
        sns.set_style("ticks")
        fig, ax = plt.subplots(1,1,figsize=(6,6)) 
   
        handles = [None] * len(reads_per_snp)
        for r in (reads_per_snp):
            cur_color = next(colors)
            cur_marker_fill_style = next(marker_fill_styles)
            cur_markeredgecolor = None
            cur_markeredgewidth = None
            if r == 2.0:
                cur_markeredgewidth = 1
                cur_markeredgecolor = 'black'
                cur_marker_fill_style = 'none'
                
            abs_error = np.sort(np.absolute(res_df.expected_dd_cfDNA[(res_df.Sim_reads_per_SNP == r) & (res_df.expected_dd_cfDNA > 0)] - res_df.OneGenome_P_donor[(res_df.Sim_reads_per_SNP == r) & (res_df.expected_dd_cfDNA > 0)])*100)
            abs_error_cdf = np.array(range(abs_error.shape[0])) * (1.0 / abs_error.shape[0])
            re_abs_error = np.sort(abs_error / (100*res_df.OneGenome_P_donor[(res_df.Sim_reads_per_SNP == r) & (res_df.expected_dd_cfDNA > 0)]) )
            #print(re_abs_error)
            ax.plot(re_abs_error, abs_error_cdf,
                    label=str(r), marker='', linestyle='-', color=cur_color, markersize=6, fillstyle=cur_marker_fill_style, markeredgecolor = cur_markeredgecolor, markeredgewidth = cur_markeredgewidth)
        #plt.xscale('log', nonposy='clip')
        #plt.yscale('log', nonposy='clip')
        plt.xlabel("Absolute relative difference between the predicted and expected precent of dd-cfDNA(\%)",fontsize=12)
        plt.ylabel("Samples with smaller or equal absolute error",fontsize=12)
        plt.legend(title="Avg. coverage",loc='lower right') 
        sns.despine()
        plt.xlim((0,5.0))
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        fig.savefig( output_plot_filename.replace("_cdf.pdf","_rel_cdf.pdf"), bbox_inches='tight')
        plt.close()
        
  

###############################################################################################################
# main function
###############################################################################################################

def main():
  
    parser = argparse.ArgumentParser("Comparing inference of donor-derived cfDNA fraction using recipient and donor genotypes and only recipient")
    
    
    parser.add_argument("org_samples_table_filename", type=str, help= "The original samples table file name")
    
    parser.add_argument("output_dir", type=str, help= "The output directory of the analysis, this directory should contain the fit directory")

    parser.add_argument("fig_output_dir", type=str, help= "the name of the directory inside the analysis directory Output/Analysis/<dir>")

    parser.add_argument("-f", "--fit_versions", dest='fit_versions', default=['V22', 'V23'], nargs='*',
                      help=('fit verions (for example: V21), can be multiple strings'))
    
    parser.add_argument("-a", "--additional_output_dirs", dest='additional_output_dirs', default=[], nargs='*',
                      help=('Additional output directories'))
    
    parser.add_argument("-c", "--compare_ver", dest='compare_ver',  action='store_true',
                        help=('compare fit versions'))
    
    
    args = parser.parse_args()
    
    output_dir = args.output_dir
    fit_versions = args.fit_versions
    fig_output_dir = output_dir + '/Analysis/' + args.fig_output_dir
    additional_output_dirs = args.additional_output_dirs

    print('Creating output directory (if does not exists)')
    
    if not os.path.exists(fig_output_dir):
        os.makedirs(fig_output_dir)

    if not os.path.exists(fig_output_dir + '/Patients'):
        os.makedirs(fig_output_dir + '/Patients')


   
    print("Drawing simulations version 0.1")
    
    samples_df_filename, res_dfs = load_simulation_results(args.org_samples_table_filename, output_dir, fit_versions, additional_output_dirs)
    
    
    
    ######################################################################
    # comparing fit results 
    ######################################################################
    
    if (args.compare_ver):
      
      for f1,fit_ver1 in enumerate(fit_versions):
        
        org_res_df = res_dfs[fit_ver1]
        
        
        # filtering samples that were not run
        res_df_all1 = org_res_df.ix[((~org_res_df['expected_dd_cfDNA'].isnull()) &
                                    (~org_res_df['OneGenome_P_donor'].isnull()) &
                                    (org_res_df.s1_SampleUniqueStr.str.contains('_BC')) &
                                    (org_res_df.s2_SampleUniqueStr.str.contains('_BC')) )  ,:].copy()
        
        for f2,fit_ver2 in enumerate(fit_versions):
          
          if (f1 < f2):
        
            org_res_df = res_dfs[fit_ver2]
            
            
            # filtering samples that were not run
            res_df_all2 = org_res_df.ix[((~org_res_df['expected_dd_cfDNA'].isnull()) &
                                        (~org_res_df['OneGenome_P_donor'].isnull()) &
                                        (org_res_df.s1_SampleUniqueStr.str.contains('_BC')) &
                                        (org_res_df.s2_SampleUniqueStr.str.contains('_BC')) )  ,:].copy()
            
            output_plot_filename = fig_output_dir + "/" + fit_ver1 + "_" + fit_ver2  + "_log" + ".pdf"
            
            plot_2_fit_ver_comparison(res_df_all1, res_df_all2, fit_ver1, fit_ver2, output_plot_filename)
            
            
      # only unrelated
      for f1,fit_ver1 in enumerate(fit_versions):
        
        org_res_df = res_dfs[fit_ver1]
        
        
        # filtering samples that were not run
        res_df_all1 = org_res_df.ix[((~org_res_df['expected_dd_cfDNA'].isnull()) &
                                    (~org_res_df['OneGenome_P_donor'].isnull()) &
                                    (org_res_df.s1_SampleUniqueStr.str.contains('_BC')) &
                                    (org_res_df.s2_SampleUniqueStr.str.contains('_BC'))  &
                                    (org_res_df['ibd_m1'] >= 50) & (org_res_df['ibd_m2'] >= 50))  ,:].copy()
        
        for f2,fit_ver2 in enumerate(fit_versions):
          
          if (f1 < f2):
        
            org_res_df = res_dfs[fit_ver2]
            
            
            # filtering samples that were not run
            res_df_all2 = org_res_df.ix[((~org_res_df['expected_dd_cfDNA'].isnull()) &
                                        (~org_res_df['OneGenome_P_donor'].isnull()) &
                                        (org_res_df.s1_SampleUniqueStr.str.contains('_BC')) &
                                        (org_res_df.s2_SampleUniqueStr.str.contains('_BC')) &
                                        (org_res_df['ibd_m1'] >= 50) & (org_res_df['ibd_m2'] >= 50))  ,:].copy()
            
            output_plot_filename = fig_output_dir + "/" + fit_ver1 + "_" + fit_ver2  + "_unrelated_log" + ".pdf"
            
            plot_2_fit_ver_comparison(res_df_all1, res_df_all2, fit_ver1, fit_ver2, output_plot_filename)
            
      exit()
            
        
        
      
    
    ######################################################################
    # predicted versus expected for non IBD for coverage 
    ######################################################################
    
    # python ~/bin/cfDNA1G/simulations/cfDNA_draw_sim_results.py ./bm.samples.tsv Output Vmulti -f V37 -a ../bm_sim15/Output
    
    for fit_ver in fit_versions:
        
        org_res_df = res_dfs[fit_ver]
        
        org_res_df.to_csv(fig_output_dir + "/" + fit_ver + ".sim.res.tsv",sep='\t',header=True,index=False)
        
        
        # filtering samples that were not run 
        res_df_all = org_res_df.ix[((~org_res_df['expected_dd_cfDNA'].isnull()) &
                                    (~org_res_df['OneGenome_P_donor'].isnull()) &
                                    (org_res_df.s1_SampleUniqueStr.str.contains('_BC')) &
                                    (org_res_df.s2_SampleUniqueStr.str.contains('_BC')) )  ,:].copy()
        
        
        res_df_all['tmp_input_conf'] = res_df_all['Sim_NumSNPs'].astype(str) + '_' + res_df_all['Sim_SNPs_selectionMethod'].astype(str) + '_' + res_df_all['Sim_reads_per_SNP'].astype(str)
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & (res_df_all['Sim_reads_per_SNP'] >  1.1) & (res_df_all['Sim_NumSNPs'] > 1000)].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_full_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_full_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        
        
        
        ##############################################################################################################################
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & ( (res_df_all['Sim_reads_per_SNP'] == 1.0) & (res_df_all['Sim_NumSNPs'] == 100) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & ( (res_df_all['Sim_reads_per_SNP'] == 1.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads1_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads1_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & ( (res_df_all['Sim_reads_per_SNP'] == 0.5) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads05_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads05_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        
        ##############################################################################################################################
       
       
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & ( (res_df_all['Sim_reads_per_SNP'] == 0.33) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads033_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsampl_reads033e_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        
        ##############################################################################################################################
       
       
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] ==2) &  (res_df_all['ibd_m2'] == 100)) & ( (res_df_all['Sim_reads_per_SNP'] < 3.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_halfsiblings_subsample_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_halfsiblings_subsample_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] ==2) &  (res_df_all['ibd_m2'] == 2)) & ( (res_df_all['Sim_reads_per_SNP'] < 3.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_siblings_subsample_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_siblings_subsample_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] ==4) &  (res_df_all['ibd_m2'] == 100)) & ( (res_df_all['Sim_reads_per_SNP'] < 3.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_firstcousins_subsample_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_firstcousins_subsample_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] ==6) &  (res_df_all['ibd_m2'] == 100)) & ( (res_df_all['Sim_reads_per_SNP'] < 3.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique()
        print("ploting related configurtions for full SNPs:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_secondcousins_subsample_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_secondcousins_subsample_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        

        ##############################################################################################################################
       
        
        
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) & ( (res_df_all['Sim_reads_per_SNP'] == 1.0) & (res_df_all['Sim_NumSNPs'] > 1000) )].unique() # (res_df_all['Sim_reads_per_SNP'] == 1.0) & 
        print("ploting related configurtions for full SNPs but partial reads:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realtness_subsample_reads_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        
        
        print("ploting unrelated subsampling cfDNA reads:")
        
        cur_res_df = res_df_all.ix[ ((res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) &  (res_df_all['Sim_NumSNPs'] > 1000) )  ,:] 
       
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_reads_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_reads_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)
        
        ##############################################################################################################################
        
        
        print("ploting unrelated subsampling cfDNA reads:")
        
        cur_res_df = res_df_all.ix[ ((res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) &  (res_df_all['Sim_NumSNPs'] ==100) )  ,:] 
       
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_reads_K100_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_reads_K100_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)
        
        
        ##############################################################################################################################
        
        
        
        
        print("ploting unrelated subsampling SNPs:")
        
        cur_res_df = res_df_all.ix[ ((res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) &  (res_df_all['Sim_reads_per_SNP'] > 1.1) )  ,:] 
       
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_SNPs_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_nonrelated_subsample_SNPs_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = True, use_markersAndFill_for_snpXread = True)
        
        
        ##############################################################################################################################
        
        
        
        tmp_input_conf = res_df_all['tmp_input_conf'][((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10))].unique()
        
        print("ploting related configurtions:")
        print(tmp_input_conf)
        
        cur_res_df = res_df_all.ix[ res_df_all.tmp_input_conf.isin(tmp_input_conf) ,:]
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realted_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realted_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realted_cdf" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "cdf", log_scale = True, only_unrelated = False)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_realted_relcdf" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "relcdf", log_scale = True, only_unrelated = False)
        
        
        print("ploting unrelated:")
        
        cur_res_df = res_df_all.ix[((res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) )  ,:] 
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_unrealted_log" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_unrealted_mae" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = True)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_unrealted_cdf" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "cdf", log_scale = True, only_unrelated = True)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "_sim_unrealted_relcdf" + ".pdf"
        plot_sim_results(cur_res_df, cur_output_plot_filename, type_of_y_axis = "relcdf", log_scale = True, only_unrelated = True)

    # old plots
    if (False):
        exit()
        
        
        cur_res_df = res_df_all.ix[ ((res_df_all.s2_SampleUniqueStr.str.contains('_BC')) & (res_df_all.s1_SampleUniqueStr.str.contains('_BC')) ),:]
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_all_log" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True, only_unrelated = False)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_all_mae" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True, only_unrelated = False)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_all_cdf" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "cdf", log_scale = True, only_unrelated = False)
        
        
        
        
        cur_res_df = res_df_all.ix[((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10) )  ,:] 
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_realted_log" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "prediction", log_scale = True)  
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_realted_mae" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "mae", log_scale = True)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_realted_cdf" + ".pdf"
        plot_sim_results_new(cur_res_df, cur_output_plot_filename, type_of_y_axis = "cdf", log_scale = True)
        
        exit()
        raise NotImplementedError("DEBUG")
        
        
        res_df = res_df_all.ix[((res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) &
                            (res_df_all['Sim_SNPs_selectionMethod'] == 'Random')  &  
                            (res_df_all['Sim_NumSNPs'] == 2000000)  )  ,:] 
        
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_all_Rand2_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results(cur_res_df, cur_output_plot_filename.replace("_log.pdf","_lin.pdf"), log_scale = False)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_all_Rand2_nonrelated_byReads_under005_log" +  ".pdf"
        cur_res_df = res_df.ix[(res_df.expected_dd_cfDNA <= 0.05) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results(cur_res_df, cur_output_plot_filename.replace("_log.pdf","_lin.pdf"), log_scale = False)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_Rand2_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results(cur_res_df, cur_output_plot_filename.replace("_log.pdf","_lin.pdf"), log_scale = False)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_Rand2_nonrelated_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results(cur_res_df, cur_output_plot_filename.replace("_log.pdf","_lin.pdf"), log_scale = False)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_W15_Rand2_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(~res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results(cur_res_df, cur_output_plot_filename.replace("_log.pdf","_lin.pdf"), log_scale = False)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_W15_Rand2_nonrelated_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (~res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        
        res_df = res_df_all.ix[( ((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) &
                            (res_df_all['Sim_SNPs_selectionMethod'] == 'Random')  &  
                            (res_df_all['Sim_NumSNPs'] == 2000000)  )  ,:] 
        
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_Rand2_related_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_Rand2_related_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
    
        #############################################################################
        
        
        res_df = res_df_all.ix[( (res_df_all['ibd_m1'] > 10) &  (res_df_all['ibd_m2'] > 10) &
                                 (res_df_all['Sim_SNPs_selectionMethod'] == 'KMostCommonInEachBlock') )   ,:]  # &   (res_df_all['Sim_NumSNPs'] == 100)  )
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_all_KMostCommon_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_all_KMostCommon_nonrelated_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[(res_df.expected_dd_cfDNA <= 0.05) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_KMostCommon_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_KMostCommon_nonrelated_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_W15_KMostCommon_nonrelated_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(~res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_W15_KMostCommon_nonrelated_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (~res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        
        res_df = res_df_all.ix[( ((res_df_all['ibd_m1'] <= 10) |  (res_df_all['ibd_m2'] <= 10)) &
                            (res_df_all['Sim_SNPs_selectionMethod'] == 'KMostCommonInEachBlock')  &  
                            (res_df_all['Sim_NumSNPs'] == 100)  )  ,:] 
        
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_KMostCommon_related_byReads_log" + ".pdf"
        cur_res_df = res_df.ix[(res_df.s2_SampleUniqueStr.str.contains('_BC')) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
                                  
        cur_output_plot_filename = fig_output_dir + "/" + fit_ver + "sim_SNPs_BC_KMostCommon_related_byReads_under005_log" + ".pdf"
        cur_res_df = res_df.ix[( (res_df.s2_SampleUniqueStr.str.contains('_BC'))  & (res_df.expected_dd_cfDNA <= 0.05) ) ,:]
        if ( cur_res_df.shape[0] > 0):
          plot_sim_results(cur_res_df, cur_output_plot_filename)
          plot_sim_results_mae(cur_res_df, cur_output_plot_filename.replace(".pdf","_mae.pdf"))
          plot_sim_results_cdf(cur_res_df, cur_output_plot_filename.replace(".pdf","_cdf.pdf"))
          
        

    
if __name__ == '__main__':
  main()

