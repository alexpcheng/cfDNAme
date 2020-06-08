#!/usr/bin/env python

#grammy-gdt -- script to parse reference genome set into GRAMMy genome data  

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

#BLD_DIR="build/lib.linux-%s-2.6/gem"
#sys.path.append( BLD_DIR % platform.machine() )
try:
  #installed import
  from grammy import gemaux, gemlib, gemcore, gemmath, gemutil
except ImportError:
  #debug import
  import gemaux, gemlib, gemcore, gemmath, gemutil
  #np.seterr(all='raise')


def main():
  #

  parser = argparse.ArgumentParser(description="grammy-gdt Commandline Tool")
  parser.add_argument("o_prefix", metavar= "o_prefix", help="output prefix, o_prefix.gdt will be the output filename")
  parser.add_argument("taxids", metavar= "taxids", help="taxids to include, taxids: t1,t2,...,tm  each tid t? is a INTEGER and must be found in grefs/gid_tid.dmp")
  parser.add_argument("-d","--dmp", dest="dmp", help="gid to tid dump file, default=grefs/gid_tid.dmp", default="grefs/gid_tid.dmp")
  parser.add_argument("-r","--ref", dest="ref", help="reference genome dir, default=grefs", default="grefs")
  parser.add_argument("-p","--per", dest="per", type=int, help="number of genomes per file, default=20", default=20)
  arg_namespace = parser.parse_args()

  taxids = [int(v) for v in vars(arg_namespace)['taxids'].split(',')]
  o_prefix = vars(arg_namespace)['o_prefix']
  dmp = vars(arg_namespace)['dmp']
  ref = vars(arg_namespace)['ref']
  per = vars(arg_namespace)['per']

  gdata = gemaux.Genome_by_gref( taxids=taxids, o_prefix=o_prefix, dmp_file=dmp, gref_dir=ref, per_set=per ) 
  gdata.write( o_prefix )

if __name__ == "__main__":
  main()
