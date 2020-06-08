#!/usr/bin/env python3
# Eilon Sharon December 2014


"""

Can get as input either single end or pair ends reads

Removing both reads from a pair if one is lower than quality threshold

"""

import sys
import getopt
import os

import argparse
import itertools

import gzip
import numpy as np


_debug = 0



def filter_reads(in_fastq_filename1,in_fastq_filename2,out_fastq_filename_pref,min_bp_q,min_read_percent):
  
  is_paired_ends =  ( len(in_fastq_filename2) > 0 and in_fastq_filename2 != "NA")
  
  if (is_paired_ends):
    print("running as paired-ends")
  else:
    print("running as single-end")
  
  fastq_in_r1 = gzip.open(in_fastq_filename1, mode='rt')
  fastq_out_r1 = gzip.open(out_fastq_filename_pref + '.R1' + '.filterred.fastq.tmp.gz', mode='wt')
  
  
  if (is_paired_ends):
    fastq_in_r2 = gzip.open(in_fastq_filename2, mode='rt')
    fastq_out_r2 = gzip.open(out_fastq_filename_pref + '.R2' + '.filterred.fastq.tmp.gz', mode='wt')

  if (is_paired_ends):
    j=0
    i=0
    for l1,l2 in zip(fastq_in_r1,fastq_in_r2):
      
      i=i+1
      
      if (i==1):
        r1=l1
        r2=l2
        
      else:
        r1 = r1 + l1
        r2 = r2 + l2
        
        if (i==4):
          i=0;
          j=j+1
          
          #print('DEBUG',type(l1),l1)
          l1  = l1.strip('\n')
          l2  = l2.strip('\n')
          
          q1 = np.array([ord(c) for c in l1])
          q2 = np.array([ord(c) for c in l2])
          
          # filter
          
          if ( sum(q1[1:]-33>min_bp_q)/(len(q1)-1) >= min_read_percent and \
               sum(q2[1:]-33>min_bp_q)/(len(q2)-1) >= min_read_percent ):
            fastq_out_r1.write(r1)
            fastq_out_r2.write(r2)
  else:
    j=0
    i=0
    for l1 in fastq_in_r1:
      
      i=i+1
      
      if (i==1):
        r1=l1
        
      else:
        r1 = r1 + l1
        
        if (i==4):
          i=0;
          j=j+1
          
          l1  = l1.strip('\n')
          
          q1 = np.array([ord(c) for c in l1])
          
          # filter
          if ( sum(q1[1:]-33>min_bp_q)/(len(q1)-1) >= min_read_percent ):
            fastq_out_r1.write(r1)
    

  fastq_in_r1.close()
  fastq_out_r1.close()

  if (is_paired_ends):
    fastq_in_r2.close()
    fastq_out_r2.close()
  
  print("renaming the tmp file to output files")
  
  os.rename(out_fastq_filename_pref + '.R1' + '.filterred.fastq.tmp.gz', out_fastq_filename_pref + '.R1' + '.filterred.fastq.gz')
  if (is_paired_ends):
    os.rename(out_fastq_filename_pref + '.R2' + '.filterred.fastq.tmp.gz', out_fastq_filename_pref + '.R2' + '.filterred.fastq.gz')
  
def usage():
  print("")
  print('Usage: ' + sys.argv[0] + ' --r1 <in fastq1> --r2 <in fastq 2> --o <output prefix> --q <min bp Q> --p <min read percent>')
  print('filtering low quality reads')
  print("")

def main():

  min_bp_q = 21
  min_read_percent = 0.50
  
  in_fastq_filename1 = ''
  in_fastq_filename2 = ''
  out_fastq_filename_pref = './tmp_'
  
  argv = sys.argv[1:]
  
  if (len(argv) == 0):
    usage()                     
    sys.exit()

  try:                                
    opts, args = getopt.getopt(argv, "h1:2:o:q:p:d", ["help","r1=","r2=","o=","q=","p=","d"])
  except getopt.GetoptError:          
    usage()                         
    sys.exit(2)                     
  for opt, arg in opts:
    if opt in ("-h", "--help"):
      usage()                     
      sys.exit()                  
    elif opt == '-d':
      global _debug               
      debug = 1                  
    elif opt in ("-1", "--r1"):
      in_fastq_filename1 = arg
    elif opt in ("-2", "--r2"):
      in_fastq_filename2 = arg
    elif opt in ("-o", "--o"):
      out_fastq_filename_pref = arg
    elif opt in ("-q", "--q"):
      min_bp_q = int(float(arg))
    elif opt in ("-p", "--p"):
      min_read_percent = float(arg)
    else:
      print('Unknown option:' + opt)
      sys.exit()


  print("Parameters:")
  print("in_fastq_filename1: " + in_fastq_filename1)
  print("in_fastq_filename2: " + in_fastq_filename2)
  print("out_fastq_filename_pref: " + out_fastq_filename_pref)
  print("min_bp_q: " + str(min_bp_q))
  print("min_read_percent: " + str(min_read_percent))
  
  
  filter_reads(in_fastq_filename1,in_fastq_filename2,out_fastq_filename_pref,min_bp_q,min_read_percent)
  
  print('Done!')


  
if __name__ == '__main__':
  main()


