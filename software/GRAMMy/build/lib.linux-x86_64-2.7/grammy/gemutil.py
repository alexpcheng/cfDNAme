import os, sys, tempfile, shutil, subprocess, time
from Bio import SeqIO
try:
  from grammy import gemlib
except ImportError:
  import gemlib

# if nprocs == 1, go this way
def simple_wet_run( rdata, gdata, bootstraps, model_par, algor ):
  # k-tuple based method, model_para = kernel(int)
  # mapping based method, model_para = mtxfile(string)
  #rdata.reads_file = '/'.join( rdata.reads_file.split('/')[1:] ) #tmp fix for t data, tbr
  fvecs = []
  llds = []
  btps = []
  us = None
  if algor == 'map':
    us = gemlib.pMONO()
    us.prep( rdata, gdata, model_par )
  if algor == 'mono':
    us = gemlib.URAE_EM()
    us.prep( rdata, gdata, model_par )
  if algor == 'dual':
    print 'model_par=', model_par
    us = gemlib.DUAL_EM()
    us.prep( rdata, gdata, model_par )
    (cdatas, C, Ci) = us.cluster()
    us.classify_all(cdatas, C, Ci)
    us.init()

  us.em.init_V('M')                                       #init by moments
  tolerance = 0.001/us.em.MN
  mle_start_time = time.time()
  us.status = us.em.solve( tolerance, 1000, 'U' )                  #in order: tolerance, max_iter, stopping rule
  mle_end_time = time.time()
  llds.append( us.em.logL )                               #log likelihood
  fvecs.append( list(us.em.V) )                                 #mixing parameter MLE
  print >>sys.stderr, "[%s] em: MLE solved in %s steps, %s seconds" % (__name__, us.status, mle_end_time - mle_start_time)
  bootstrap_start_time = time.time()
  btps.append( us.em.bootstrap( bootstraps, tolerance, 1000, 'M', 'L' ) )#in order: bootstrap_times, tolerance, max_iter, init_method, stopping rule
  bootstrap_end_time = time.time()
  print >>sys.stderr, "[%s] em: Bootstrapped %s times, %s seconds" % (__name__, bootstraps, bootstrap_end_time - bootstrap_start_time)
  return fvecs, llds, btps, us

def fasta_reads_file( reads_file ):
  if not os.path.exists( reads_file ):
    print >>sys.stderr, ' reads_file not found '
    exit()
  base_reads_file = os.path.basename( reads_file )
  tmp_reads_file = reads_file  # set to be original file

  if base_reads_file.endswith('.gz'):
    #decompress reads_file
    TMP_DIR = tempfile.mkdtemp()
    tmp_reads_file = os.path.join(TMP_DIR, base_reads_file)
    shutil.copyfile( reads_file, tmp_reads_file )
    args = 'gzip -d ' + tmp_reads_file
    subprocess.Popen( args, shell=True).communicate()
    tmp_reads_file = os.path.splitext( tmp_reads_file )[0] # set to be decompressed file

  #print os.path.splitext( tmp_reads_file )[1]
  if os.path.splitext( tmp_reads_file )[1] in [ '.fq', '.fastq' ]:
    #convert to fasta file
    old_tmp_reads_file = tmp_reads_file
    tmp_reads_file = os.path.splitext( old_tmp_reads_file )[0] + '.fasta' # set to be converted file
    SeqIO.convert( old_tmp_reads_file, 'fastq', tmp_reads_file, 'fasta' )

  return tmp_reads_file

def blast8_hits_condense( hits_set, rl ): #simple greedy condense
  merged_set = []
  rest_hits = hits_set
  while rest_hits != []:
    merged_hit, rest_hits = merge_hits( rest_hits, rl )
    #print merged_hit, rest_hits
    merged_set.append( merged_hit )
  return merged_set

def merge_hits( hits_set, rl ):
  basis_hit = hits_set.pop(0)
  to_remove = []
  for hit in hits_set:
    if hit_compatible( basis_hit, hit, rl ):
      basis_hit = hit_merge( basis_hit, hit )
      to_remove.append( hit )
  for hit in to_remove:
    hits_set.remove( hit )
  return basis_hit, hits_set

def hit_compatible( basis_hit, hit, rl ):
  # check if hits are compatible
  if basis_hit[-1] != hit[-1]:
    #print "not the same direction!"
    return False
  else:
    #print 
    rseg = hit[0][0]
    gseg = hit[1][0]
    for seg in basis_hit[0]:
      #print basis_hit[0], 
      if not ( seg[1] < rseg[0] or seg[0] > rseg[1] ): #
        #print "read overlapping!"
        return False
    for seg in basis_hit[1]:
      if not ( ( ( seg[1] < gseg[0] or seg[0] > gseg[1] ) and basis_hit[-1] ) or ( ( seg[1] > gseg[0] or seg[0] < gseg[1] ) and ( not basis_hit[-1] ) ) ):
        #print "genome overlapping!"
        return False
    gmin = 100000000000
    gmax = 0
    for seg in basis_hit[1]:
      for tip in seg:
        if tip < gmin:
          gmin = tip
        if tip > gmax:
          gmax = tip
    if not ( ( basis_hit[-1] and ( gmax-rl <= gseg[0] and gmin+rl >= gseg[1] ) ) or ( ( not basis_hit[-1] ) and ( gmax-rl <= gseg[1] and gmin+rl >= gseg[0] ) ) ):
      #print "genome out of range!"
      return False
  return True

def hit_merge( basis_hit, hit ):
  #      (rs,           re)               (gs,             ge)              id               al               ev                direction: T=same, F=reverse
  basis_hit[0].append( hit[0][0] )
  basis_hit[1].append( hit[1][0] )
  basis_hit[3] = basis_hit[3] + hit[3]
  if hit[4] < basis_hit[4]:
    basis_hit[4] = hit[4]
  basis_al = 0
  for seg in basis_hit[0]:
    basis_al += abs( seg[0] - seg[1] )
  basis_hit[2] = ( (basis_hit[2] * basis_al) + hit[2] * abs(hit[0][0][1]-hit[0][0][0]) + 1 ) / ( basis_al + abs(hit[0][0][1] - hit[0][0][0]) + 1)
  return basis_hit
