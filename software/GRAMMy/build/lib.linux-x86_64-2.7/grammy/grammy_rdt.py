#!/usr/bin/env python
#grammy-rdt -- script to parse read set into GRAMMy read data  

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

import os, sys, platform, shutil, gzip, argparse
from Bio import SeqIO

#BLD_DIR="build/lib.linux-%s-2.6/gem"
#sys.path.append( BLD_DIR % platform.machine() )

try:
  #installed import
  from grammy import gemaux, gemlib, gemcore, gemmath, gemutil
except ImportError:
  #debug import
  import gemaux, gemlib, gemcore, gemmath, gemutil
  #np.seterr(all='raise')


"""suppose we have an o_prefix dir,
   all file paths are relative from that,
   rdt file is to be kept in o_prefix dir"""

def parse_reads( o_prefix, src_dir, i_suffix, name_change, seq_tech ):
  #
  for file in os.listdir( os.path.join( src_dir ) ):
    #
    if file.endswith(i_suffix):
      #
      print >>sys.stderr, "processing", file
      reads_src = None
      set_name = os.path.basename(file).rstrip(i_suffix)
      set_name = set_name in name_change and name_change[set_name] or set_name
      reads_src = os.path.join( src_dir, file )
      reads_tgt = os.path.join( o_prefix, '%s.fasta.gz' % set_name )
      if file.endswith('.gz'):
        #
        #print reads_src, reads_tgt
        shutil.copy2( reads_src, reads_tgt )
      else:
        #
        reads_content = open( reads_src, 'rb' )
        reads_tgt_file = gzip.open( reads_tgt, 'wb' )
        reads_tgt_file.writelines( reads_content )
        reads_tgt_file.close()
      reads_len = 0
      reads_no = 0
      for seq_rec in SeqIO.parse( gzip.open( reads_tgt ), "fasta" ):
        #
        reads_no += 1
        reads_len += len(seq_rec.seq)
      reads_file = '%s.fasta.gz' % set_name
      rdata = gemaux.Read_Data()
      rdata.read_tech = seq_tech
      rdata.read_length = reads_len/reads_no
      rdata.reads_file = reads_file
      rdata.reads_number = reads_no
      rdata.write( os.path.join(o_prefix,set_name) )

def main():
  #
  parser = argparse.ArgumentParser(description="grammy-rdt Commandline Tool")

  parser.add_argument("i_prefix", metavar= "i_prefix", help="itput dir prefix, a dir where reads files reside")
  parser.add_argument("o_prefix", metavar= "o_prefix", help="output dir prefix, the output will be o_prefix/xxx.rdt, use '.' for current dir")
  parser.add_argument("-s","--suf", dest="suf", help="read files suffix, default=fa.gz", default="fa.gz")
  parser.add_argument("-t","--tec", dest="tec", help="sequencing tech, default=sanger", default="sanger")
  parser.add_argument("-c","--chg", dest="chg", help="name change set 'o1:n1,o2:n2', default= ", default="")
  arg_namespace = parser.parse_args()

  i_prefix = vars(arg_namespace)['i_prefix']
  o_prefix = vars(arg_namespace)['o_prefix']
  chg = vars(arg_namespace)['chg']
  tec = vars(arg_namespace)['tec']
  suf = vars(arg_namespace)['suf']

  if chg:
    name_change = dict( [ v.split(':') for v in chg.split(',') ] )
  else:
    name_change = dict()

  parse_reads( o_prefix=o_prefix, src_dir=i_prefix, i_suffix=suf, name_change=name_change, seq_tech=tec )

if __name__ == "__main__":
  #
  main()
