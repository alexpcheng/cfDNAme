#!/usr/bin/env python

#grammy-post -- script to parse GRAMMy EM calculations after grammy-em 

#License: BSD

#Copyright (c) 2008 Li Charlie Xia
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

import os, sys, platform, argparse
import numpy as np
from scipy.io import mmio
from scipy import sparse

#BLD_DIR="./build/lib.linux-%s-2.6/gem"
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

  parser = argparse.ArgumentParser(description="grammy-post Commandline Tool")
  parser.add_argument("mix_par", metavar= "mix_par", help="untransformed mixing parameter estimates, from grammy-em")
  parser.add_argument("gen_dat", metavar= "gen_dat", help="gen_dat.gdt genome data file, from grammy-gdt")
  parser.add_argument("btp", metavar= "btp", help="bootstrap file, from grammy-em")
  arg_namespace = parser.parse_args()

  mix_par = vars(arg_namespace)['mix_par']
  gen_dat = vars(arg_namespace)['gen_dat']
  btp = vars(arg_namespace)['btp']

  gdt = gemaux.Genome_Data()
  gdt.read( gen_dat )
  f = np.loadtxt( mix_par )
  bts = np.loadtxt( btp )
  l = np.array( [ gdt.taxid_length[tid] for tid in gdt.taxids ], dtype='float' )
  #print l, f
  assert len(l) == len(f) - 1
  a = gemmath.inverse_proportion_normalize( f, l ) 
  abd_bsd = gemmath.bootstrap_standard_error( f, l, bts )
  est_avl = gemmath.weighted_average( a, l )

  #output >> .gra
  #output >> .avl
  prj,suffix = os.path.splitext(mix_par)
  gra_file = open( ".".join([prj,"gra"]), 'w' )
  print >>gra_file, '\n'.join( [ '\t'.join( [ str(v) for v in gdt.taxids ] ), '\t'.join( [ "%.4g" % v for v in a ] ), '\t'.join( [ "%.4g" % v for v in abd_bsd ] ) ] )  
  avl_file = open( ".".join([prj,"avl"]), 'w' )
  print >>avl_file, est_avl

if __name__ == "__main__":
  #
  main()
