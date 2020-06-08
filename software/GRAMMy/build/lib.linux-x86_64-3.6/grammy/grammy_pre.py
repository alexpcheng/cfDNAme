#!/usr/bin/env python

#grammy-pre -- preparing for GRAMMy EM calculation 

#License: BSD

#Copyright (c) 2008 Li Charles Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os, sys, tempfile, platform, gzip, argparse
from Bio.Blast import NCBIStandalone
from Bio import SeqIO
from scipy.io import mmio
from scipy import sparse
import scipy as sp

TMP_DIR = tempfile.mkdtemp()
#BLD_DIR="build/lib.linux-%s-2.6/gem"
#sys.path.append( BLD_DIR % platform.machine() )
#from grammy import gemaux, gemlib, gemcore, gemmath, gemutil
try:
  #installed import
  from grammy import gemaux, gemlib, gemcore, gemmath, gemutil
except ImportError:
  #debug import
  import gemaux, gemlib, gemcore, gemmath, gemutil
  #np.seterr(all='raise')

def main():
  #
#  usage_error = """ type %s -h for more information """
#  usage = """usage: %prog [options] read_dat gen_dat
#
#Arguments:
#  read_dat  STRING        will use read_dat.rdt as read data file
#  gen_dat   STRING        will use gen_dat.gdt as genome data file"""

  parser = argparse.ArgumentParser(description="grammy-pre Commandline Tool")
  parser.add_argument("read_dat", metavar= "read_dat", help="will use read_dat.rdt as read data file")
  parser.add_argument("gen_dat", metavar= "gen_dat", help="will use gen_dat.gdt as genome data file")
  parser.add_argument("-m","--mtd", dest="mtd", choices=['tbl','bam'],
    help="method for read assignment,\n  tbl -- tabular blast format",
    default='tbl')
  parser.add_argument("-p","--par1", dest="par1",
    help="first parameter for read assignment method, tbl filename or k",
    default='read_dat.tblat')
  parser.add_argument("-q","--par2", dest="par2", 
    help="second parameter for read assignment method, al,id,ev or bg",
    default="75,75,-5") 
  arg_namespace = parser.parse_args()

  rdt_prefix = os.path.basename(vars(arg_namespace)['read_dat'])
  rdt = gemaux.Read_Data()
  rdt.read( rdt_prefix )

#  gdt_prefix = os.path.basename(vars(arg_namespace)['gen_dat'])
  gdt_prefix = vars(arg_namespace)['gen_dat']
  gdt = gemaux.Genome_Data()  
  gdt.read( gdt_prefix )

  mtd = vars(arg_namespace)['mtd']
  par1 = vars(arg_namespace)['par1']
  par2 = vars(arg_namespace)['par2']
  
  if mtd == 'tbl':

    #construct dict{ idx: id }
    read_genomes = dict()
    rid_ri = dict()
    rid_rl = dict()
    i = 0

    for seq in SeqIO.parse( gzip.open( rdt.reads_file ), 'fasta' ):
      rid_ri[ seq.id ] = i
      rid_rl[ seq.id ] = len(seq.seq)
      i += 1
    reads_number = i
    kept = 0

  
    #indegrients
    if par1 == 'read_dat.tblat':
      tbl_filename = '%s.tblat' % rdt_prefix
    else:
      tbl_filename = par1
    tbl_pars = [ int(v) for v in par2.split(',') ]
    mtx_file = rdt_prefix+".mtx"

    #parsing parameters
    ALG_LEN_THRESH = tbl_pars[0]
    ALG_ID_THRESH = tbl_pars[1]
    E_VALUE_THRESH = 0.1**tbl_pars[2]

    #parsing
    hits_set = []
    for i in range(1,len(gdt.genomes_file)+1):
      #
      tbl_file = (tbl_filename+".%d") % i
      current_gi = None
      current_qi = None
      map_handle = open( tbl_file, 'rU' )
      for line in map_handle:
        #
        line = line.rstrip('\n')
        cells = line.split('\t')
        qi = rid_ri[ cells[0] ]
        gi = int(cells[1].split('|')[1])
        if current_gi != gi or current_qi != qi:
          #
          if hits_set != []:
            #print >>sys.stderr, hits_set
            #print >>sys.stderr, "pre --------"
            hits_set = gemutil.blast8_hits_condense( hits_set, rid_rl[cells[0]] )
            #print >>sys.stderr, hits_set
            #print >>sys.stderr, "post -----------------------------------------"
            for hit in hits_set:
              if hit[4] < E_VALUE_THRESH and hit[3] >= ALG_LEN_THRESH and hit[2] >= ALG_ID_THRESH:
                  #retain hits pass a E-value, total aligned length and identity Threshold
                  #print "qi=", qi, ",gi=", gi, ",gdt.gi_taxid[gi]=",gdt.gi_taxid[gi],",tid_index=", \
                  #      gdt.taxids.index(gdt.gi_taxid[gi])  
                read_genomes.setdefault( qi, [] ).append( gdt.taxids.index(gdt.gi_taxid[gi]) ) #adjust for -1 difference
                kept += 1
          hits_set = []
          current_gi = gi
          current_qi = qi
        hits_set.append(
          [ [(int(cells[6]), int(cells[7]))], [(int(cells[8]), int(cells[9]))], float(cells[2]), float(cells[3]), float(cells[10]), int(cells[8])<int(cells[9]) ]
        )
        #      rs           re               gs             ge              id               al               ev                direction: T=same, F=reverse
        #print >>sys.stderr, "Number of Total Hits=", kept
    smm = sparse.lil_matrix( (reads_number, len(gdt.taxids)+1) ) # numpy matrix is addressed from 0, we need allocate one more for unknown
    tnm = 0
  
    #first filling in mapped reads
    mn = len(read_genomes)
    for r in read_genomes:
      #
      counts = {}
      nm = len(read_genomes[r])
      for g in read_genomes[r]:
        #
        counts[g] = counts.get(g,0)+1 # fast counts unique elements in a list
      for g in counts:
        #smm[r,g] = float( counts[g] )/nm  # simplest counts and divide probability, perhaps E-value based, weighted probability 
        #smm[r,g] = float( counts[g] )/nm  # simplest counts and divide probability, perhaps E-value based, weighted probability 
        #print float( counts[g] )
        #print gdt.taxid_length[ gdt.taxids[g] ]  # simplest counts and divide probability, perhaps E-value based, weighted probability 
        smm[r,g] = float( counts[g] )/gdt.taxid_length[ gdt.taxids[g] ]  # simplest counts and divide probability, perhaps E-value based, weighted probability 
      tnm += nm
    assert tnm==kept
  
    #second imputing in unmapped reads with the average of mapped reads
    unknown_impute = smm.sum()/smm.nnz
    unknown_index = len( gdt.taxids )  # index starts from 0
    unmapped_reads = 0
    for r in range(0, reads_number):
      #
      if r not in read_genomes:
        #
        smm[r,unknown_index] = unknown_impute
        unmapped_reads += 1
    print >>sys.stderr, "generating matrix market file(.mtx)..."
    sp.io.mmwrite( mtx_file, smm, comment='E_VALUE_THRESH=%s, ALG_LEN_THRESH=%s, MAPPED=%s, A_MAPPED=%s' \
        % (str(E_VALUE_THRESH), str(ALG_LEN_THRESH), str(mn), str(smm.nnz-mn)), field='real', precision=8 )
  
  elif mtd == "kmer":
    #
    kmer_k = int(par1)
    kmer_bg = par2
    print >>sys.stderr, "Sorry, not integrated to software release yet!"


if __name__ == "__main__":
  #
  main()
