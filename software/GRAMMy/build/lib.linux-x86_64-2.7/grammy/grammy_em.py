#!/usr/bin/env python

#grammy-em -- script to carry out GRAMMy EM calculation 

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


import os, sys, platform, time, argparse

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
#  usage_error = """ type %s -h for more information """ % os.path.basename(sys.argv[0])
#  usage = """ %prog [options] read_prob_file.mt1,read_prob_file.mt2,...
#Arguments:
#  read_prob_file: input probability matrix file"""

  parser = argparse.ArgumentParser(description="grammy-em Commandline Tool")
  parser.add_argument("read_prob_file", metavar= "read_prob_file", help="input read probability matrix file from grammy-pre")
  parser.add_argument("-b","--btp", dest="btp", type=int, help="bootstrap number, default=10", default=10)
  parser.add_argument("-t","--tol", dest="tol", type=float, help="tolerance for stopping, default=10e-6", default=.000001)
  parser.add_argument("-c","--mtd", dest="mtd", help="convergenece method, (U)niform, (L)ikelihood, default=U", choices=['U','L'], default='U')
  parser.add_argument("-n","--mit", dest="mit", type=int, help="maximum number of iteration, default=1000", default=1000)
  parser.add_argument("-i","--ini", dest="ini", help="initilization method, (M)oment, (R)andom, default=M", choices=['M','R'], default='M')
  arg_namespace = parser.parse_args()

  read_prob_file = vars(arg_namespace)['read_prob_file']
  btp = vars(arg_namespace)['btp']
  tol = vars(arg_namespace)['tol']
  mtd = vars(arg_namespace)['mtd']
  mit = vars(arg_namespace)['mit']
  ini = vars(arg_namespace)['ini']

  gem = gemlib.pMONO()

  #mtx, btp, tol, mtd, mit, ini
  gem.prep( rdt=None, gdt=None, mmf=read_prob_file.split(',') )                      #mmf is a list of matrix market files to be read in
  gem.em.init_V(ini)                                                          #init by moments
  tol = tol/gem.em.MN
  mle_start_time = time.time()                                                #time mark for MLE
  gem.status = gem.em.solve( tol, mit, mtd )                                  #in order: tolerance, max_iter, stopping rule
  mle_end_time = time.time()
  lld = gem.em.logL                                                           #log likelihood
  fvec = list(gem.em.V)                                                        #mixing parameter MLE
  print >>sys.stderr, "[%s] em: MLE solved in %s steps, %s seconds" % (sys.argv[0], gem.status, mle_end_time - mle_start_time)
  bootstrap_start_time = time.time()
  btps = gem.em.bootstrap( btp, tol, mit, ini, mtd )
  #in order: bootstrap_times, tolerance, max_iter, init_method, stopping rule
  bootstrap_end_time = time.time()
  print >>sys.stderr, "[%s] em: Bootstrapped %s times, %s seconds" % (sys.argv[0], btp, bootstrap_end_time - bootstrap_start_time)

  #output 
  #print os.path.splitext(args[0])
  prj,suffix = os.path.splitext(read_prob_file)
  btp_fn = '.'.join( [prj, "btp"] )
  est_fn = '.'.join( [prj, "est"] )
  lld_fn = '.'.join( [prj, "lld"] )
  print >>open(btp_fn, 'w'), '\n'.join( [ '\t'.join(  [ "%.4g" % v for v in b ] ) for b in btps ] ) 
  print >>open(lld_fn, 'w'), lld
  print >>open(est_fn, 'w'), '\t'.join( [ "%.4g" % v for v in fvec ] )

if __name__ == "__main__":
  main()
