#!/usr/bin/env python3
# Eilon Sharon Oct 2014


"""
test whether Illumina SNP match the orientation illumina claim it is

"""

#from __future__ import division


import sys
import os
import pysam
import argparse
from Bio import SeqIO
import pandas as pd
import re
	
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein




def main():

  if len(sys.argv) != 4:
    print('usage: filter_samfile_2_bamfile.py illumina_manifest ref_fa_file output_filenme file_contains_intensity')
    sys.exit(1)

  
  illumina_manifest = sys.argv[1]
  ref_fa_file = sys.argv[2]
  output_filenme = sys.argv[3]
  
  f_out = open(output_filenme,'w')
  f_out.write("IlmnID\tName\tSNP\tAllele1\tAllele2\tChr\tPosition\tright_strand\tis_top\tShouldFlip\tcur_topbot\ttop_strand\tis_flip_top\n")

  # column data types
  
  #col_datatypes = {'IlmnID' : str, 'Name' : str,
  #                 'IlmnStrand' : str,'SNP' : str,
  #                 'AddressA_ID' : str, 'AlleleA_ProbeSeq' : str,
  #                 'AddressB_ID' : str, 'AlleleB_ProbeSeq' : str,
  #                 'GenomeBuild' : str, 'Chr' : str, 'MapInfo' : str,
  #                 'Ploidy' : str, 'Species' : str, 'Source' : str,
  #                 'SourceVersion' : str, 'SourceStrand' : str,
  #                 'SourceSeq' : str, 'TopGenomicSeq' : str, 'BeadSetID' : str,
  #                 'Exp_Clusters' : str, 'Intensity_Only' : str, 'RefStrand' : str}
  
  sys.stderr.write('Loading illumina manifest: %s\n' % (illumina_manifest))
  illu_man_df =  pd.read_csv(illumina_manifest, skiprows = 7, dtype = {'AlleleB_ProbeSeq' : str, 'Chr' : str })
  
  cotrols_lst = illu_man_df[illu_man_df['IlmnID'] == '[Controls]'].index.tolist()
  if len(cotrols_lst)  > 0:
    illu_man_df = illu_man_df[0:cotrols_lst[0]]
  
  #with open(illumina_manifest) as f:
  #  rows = f.readlines()
    
  #sys.stderr.write('Num rows: %d\n' % (len(rows)))
  sys.stderr.write('Num rows: %d\n' % (illu_man_df.shape[0]))

  handle = open(ref_fa_file, "rU")
  record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
  handle.close()
  
  all_cnt = 0;
  error_cnt = 0;
  wrong_strand_cnt = 0
  ok_cnt = 0
  bad_align_cnt = 0
  
  ref_strand_col_i = 20
  
  for idx,row in illu_man_df.iterrows():
      
    #sys.stderr.write('parsing line: %d\n' % (idx))
    if ((idx % 1000) == 0):
      sys.stderr.write('Parsing line: %d\n' % (idx))


    all_cnt = all_cnt+1
    
    cur_snp = row['SNP']
    cur_topbot = row['IlmnStrand']
    cur_chr_name =  'chr' + row['Chr']
    cur_chr_coor = int(float(row['MapInfo']))
    cur_starnd = row['RefStrand']
    cur_TopGenomicSeq_org = row['TopGenomicSeq']
    cur_SourceSeq = row['SourceSeq']
    
    if (cur_chr_name not in record_dict):
      error_cnt=error_cnt+1
      sys.stderr.write("ERROR1 unknow chr: %s\n" % (cur_chr_name))
      sys.stderr.write("Row: %s\n\n" % (str(row)))
      continue
    
    if (cur_chr_coor < 1):
      error_cnt=error_cnt+1
      sys.stderr.write("ERROR2 zero coordinate: %d\n\n" % (cur_chr_coor))
      continue
    
    
    cur_TopGenomicSeq = cur_TopGenomicSeq_org.upper()
    
    m = re.match(r"([a-zA-Z]+)\[(\w)/(\w)\]([a-zA-Z]+)",cur_TopGenomicSeq)

    if (not m):
      error_cnt=error_cnt+1
      continue
  
    #cur_seq = m.group(1) + m.group(2) + m.group(4)
    cur_seq = m.group(1) + 'N' + m.group(4)
  
    cur_start_coor = cur_chr_coor - 1 - len(m.group(1)) - 1 + 1
    cur_end_coor = cur_start_coor + len(cur_seq)
  
    
    cur_chr = record_dict[cur_chr_name]
  
    cur_chr_seq = cur_chr.seq[cur_start_coor:cur_end_coor]
    cur_chr_seq = cur_chr_seq.upper()
  
    cur_chr_seq_str = str(cur_chr_seq)
    cur_chr_seq_str_l = list(cur_chr_seq_str)
    cur_chr_seq_str_l[len(m.group(1))] = 'N'
    cur_chr_seq_str = ''.join(cur_chr_seq_str_l)
    
    # switching to N any non standard
    alpha = ['A','T','C','G']
    
    cur_chr_seq_str2 = ''.join([i if (i in alpha and j in alpha) else 'N' for i,j in zip(cur_chr_seq_str,cur_seq)])
    cur_seq2         = ''.join([j if (i in alpha and j in alpha) else 'N' for i,j in zip(cur_chr_seq_str,cur_seq)])
  
  
  
    #cur_chr_seq = Seq(cur_chr_seq_str2, generic_dna)
    
    cur_chr_dna = Seq(cur_chr_seq_str2, generic_dna) # Seq(cur_chr_seq, generic_dna)
    cur_seq_dna = Seq(cur_seq2, generic_dna)

    is_plus  = str(cur_chr_dna) == str(cur_seq_dna)
    is_minus = str(cur_chr_dna) == str(cur_seq_dna.reverse_complement())
  
    if (cur_starnd == '+'):
      is_ref_plus = True
    else:
      is_ref_plus = False
  
    
    if (cur_topbot == 'TOP'):
      is_top = True
    else:
      is_top = False
      
      
    
    
    m = re.match(r"\[(\w)/(\w)\]",cur_snp)
    
    if (not m):
      error_cnt=error_cnt+1
      #sys.stderr.write("%d ERROR3 string (not SNP):\n%s\n%s\n\n" % (error_cnt, cur_TopGenomicSeq, row))
      continue
    
    cur_allele1 = m.group(1)
    cur_allele2 = m.group(2)

  
    if ( (is_ref_plus and is_plus and is_top) or \
         (is_ref_plus and is_minus and not is_top) or \
         (not is_ref_plus and is_minus and is_top) or \
         (not is_ref_plus and is_plus and not is_top) ):
      is_ok = True
      is_wrong_strand = False
    elif ( (not is_ref_plus and is_plus and is_top) or \
           (not is_ref_plus and is_minus and not is_top) or \
           (is_ref_plus and is_minus and is_top) or \
           (is_ref_plus and is_plus and not is_top) ):
      is_ok = True
      is_wrong_strand = True
      wrong_strand_cnt=wrong_strand_cnt+1
    else:
      is_ok = False
      is_wrong_strand = False

    if (is_ok):
      ok_cnt=ok_cnt+1
    else:
      bad_align_cnt=bad_align_cnt+1
  
    if (not is_ok or is_wrong_strand): #not is_ok or 
      #cur_chr = 
      sys.stderr.write("%d , %d ____________________________________________________________________\n" % (idx, wrong_strand_cnt))
      sys.stderr.write('++++++++++' + row['SNP'] + ' alleles ' + cur_allele1 + ',' + cur_allele2 + '\n')
      sys.stderr.write('++++++++++' + row['IlmnStrand'] + '\n')
      sys.stderr.write('++++++++++' + cur_chr_name + '\n')
      sys.stderr.write('++++++++++' + ("%d" % cur_chr_coor) + '\n')
      sys.stderr.write('++++++++++:::' + cur_starnd + '\n')
      #sys.stderr.write('++++++++++ source:' + cur_SourceSeq + '\n')
      sys.stderr.write('++++++++++    TOP:' + cur_TopGenomicSeq_org + '\n')
      sys.stderr.write('++++++++++cur seq:' + str(cur_seq_dna) + '\n')
      sys.stderr.write('++++++++++seq chr:' + str(cur_chr_dna) + '\n')
      sys.stderr.write('++++++++++\n')
      sys.stderr.write('++++++++++cur seq:' + str(cur_seq_dna) + '\n')
      sys.stderr.write('++++++++++rev chr:' + str(cur_chr_dna.reverse_complement()) + '\n')
  
      sys.stderr.write('++++++++++is   ref plus:' + str(is_ref_plus) + '\n')
      sys.stderr.write('++++++++++is match plus:' + str(is_plus) + '\n')
      sys.stderr.write('++++++++++is match mins:' + str(is_minus) + '\n')
      sys.stderr.write('++++++++++is        top:' + str(is_top) + '\n')
      sys.stderr.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' + '\n')
      sys.stderr.write('++++++++++is bad strand:' + str(is_wrong_strand) + '\n')
      sys.stderr.write('++++++++++is         ok:' + str(is_ok) + '\n')
    
      sys.stderr.write('____________________________________________________________________' + '\n' + '\n')
    
    
    
    right_strand = cur_starnd
    
    if (is_wrong_strand):
      if (cur_starnd == '-'):
        right_strand = '+'
      else:
        right_strand = '-'
        
    top_strand = right_strand
    if (not is_top):
      if (top_strand == '-'):
        top_strand = '+'
      else:
        top_strand = '-'
        
    is_flip_top = False
    if (top_strand == '-'):
      is_flip_top = True
    
    #f_out.write("IlmnID\tName\tSNP\tAllele1\tAllele2\tChr\tPosition\tright_strand\tis_top\tShouldFlip\tcur_topbot\ttop_strand\tis_flip_top\n")

    
    if (is_ok):
      f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\n" % \
                  (row['IlmnID'],row['Name'],row['SNP'], \
                   cur_allele1,cur_allele2,cur_chr_name,\
                   cur_chr_coor, right_strand,is_top,(right_strand=='-'), \
                   cur_topbot, top_strand, is_flip_top) )
  
  sys.stderr.write('____________________________________________________________________\n')
  sys.stderr.write('____________________________________________________________________\n')
  sys.stderr.write('cnt          all: %d\n' % (all_cnt))
  sys.stderr.write('cnt        error: %d\n' % (error_cnt))
  sys.stderr.write('cnt wrong strand: %d\n' % (wrong_strand_cnt))
  sys.stderr.write('cnt           OK: %d\n' % (ok_cnt))
  sys.stderr.write('cnt    bad align: %d\n' % (bad_align_cnt))
  sys.stderr.write('____________________________________________________________________\n')
  sys.stderr.write('____________________________________________________________________\n\n')
    
  f_out.close()
  
  
  

if __name__ == '__main__':
  main()
 