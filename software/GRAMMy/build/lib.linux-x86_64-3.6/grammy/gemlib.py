import os, sys, subprocess, time, math, random
import numpy as np
import scipy as sp
import scipy.io
from Bio import SeqIO

try:
  from grammy import gemcore
except ImportError:
  import gemcore

def gclassify_one(cdata, C, gdt, genomes_records, L, K, i):

  # classification function for DUAL_EM

  # mpi controls are stripped away, left to control script
  mc = gemcore.VectorFloat( C, 0 )

  #print "[STATUS: classifying gi=%d" % i
  start_time = time.time()
  gemcore.gclassify(C, cdata, \
    'N'.join( [ genomes_records[gi].seq.tostring() for gi in gdt.taxid_genomes[gdt.taxids[i]] ] ), \
        L, K, mc)
  end_time = time.time()
  print "one classification, total time: %d" % (end_time - start_time)
    
  return mc

def sample_wr(population, k):
  "Chooses k random elements (with replacement) from a population"
  n = len(population)
  _random, _int = random.random, int  # speed hack 
  result = [None] * k
  for i in xrange(k):
    j = _int(_random() * n)
    result[i] = population[j]
  return result

def rand_f(k):
  f = np.random.random(k)
  return f/f.sum()

class pMEM(object): #python base class for cMEM engine
  em = None
  M = None
  N = None
  gl = None
  rl = None
  s = None
  f = None
  a = None
  status = None
  egl = None
  gdt = None

  def seed(self):
    self.s = rand_f(self.M)
    self.em.V = gemcore.VectorFloat(self.s)
    #print >>sys.stderr, 'seeded, new V:', list(self.em.V)

  def reset(self, R=100.):
    self.em.R = R   # raise residue to run again
    self.em.del_P() # hash matrix need to be renewed
    self.em.get_P()
    print >>sys.stderr, "reset: %f => R" % self.em.R

  def solf(self, c=0.01, n=1000, c_mode='L' ): # 0.0001 for 1% of 1%, 0.01 for efficient solver
    c = 1./float(self.M)*c      
    # automatic convergen error criterion detection
    # convergence error would affect relative error at most by 1% of each percent of the relative error 
    self.status = self.em.solve( c, n, c_mode )
    self.f = np.array( list(self.em.V) )
    #print "f=", '\t'.join([ '%.5g' % f for f in self.f ])
    #print "egl=", '\t'.join([ '%s' % l for l in self.egl ])

  def sola(self):
    M = self.M
    f = self.f
    egl = self.egl
    c = 0
    for j in range(0, M):
      c += f[j]/egl[j]
    #print "c=", c
    self.a = c*f/egl
    self.a = self.a/sum(self.a)
    #print "a=", '\t'.join([ '%.5g' % a for a in self.a ])

  #def em_free(self):
  #  gemcore.delocate( self.em.P, self.N )

  #def prep
  #def prep_sim
  #above to be overrided

class pMONO(pMEM):
  #Mapping Based Estimation
  mmf = []
  mmf_info = {}
  #egl = []

  def get_mmf_info(self): #require self.mmf already set and in list
    self.N = None
    self.M = None
    for file in self.mmf:
      try:
        tmp_info = sp.io.mmio.mminfo( file )
      except IOError:
        sys.exit("exit from pMONO, can't open mmf file: %s!" % file)
      self.mmf_info[file] = [ int(v) for v in tmp_info[:3] ] # N, M, NNZ in order
      if not self.M:
        self.M = self.mmf_info[file][1] 
      else:
        assert self.M == self.mmf_info[file][1]
      if not self.N:
        self.N = self.mmf_info[file][0]
      else:
        self.N += self.mmf_info[file][0]

  def init_Prg( self ):
    offset = 0
    for file in self.mmf:
      self.em.init_Prg(offset, file)
      offset += self.mmf_info[file][0]

  # prepare in working mode, run seed() first
  def prep(self, rdt=None, gdt=None, mmf=None): 
    """ rdt is an object of Read_Data type, 
        required: rdt.reads_files, rdt.read_length
        gdt is an object of Genome_Data type
        required: gdt.genomes_file, gdt.taxid_genomes
        mmf is a list of mapping matrix of scipy sparse matrix for matrix type
        in the order as reads_files
        required: normalized, each mm(i,j) = Prob( ri \in gj ) """
    
    ### Dealing mmf files
    self.mmf = mmf
    self.get_mmf_info()

    """
    ### Dealing gdt, outdated
    print >>sys.stderr, "using genome file:", gdt.genomes_file

    if os.environ.get("TMPDIR") != None:
      print >>sys.stderr, "solving concurrent gdt.genomes_file access problem by using %s on this node" % os.environ.get("TMPDIR") 
      r = subprocess.Popen( "cp " + gdt.genomes_file + " $TMPDIR/genomes.fna", shell=True ).communicate()
      if r != (None, None):
        print >>sys.stderr, "return =", r, ", error!"
        print >>sys.stderr, "Error copying temporary genomes.fna file, quitting..."
        sys.exit("exit form pMONO...")
      gdt.genomes_file = os.environ.get("TMPDIR") + "/genomes.fna"
    print >>sys.stderr, "proceed to read genomes.fna file"

    genomes_handle = open( gdt.genomes_file, "r" )
    genomes_records = []
    for seq_record in SeqIO.parse(genomes_handle, "fasta"):
      genomes_records.append(seq_record)

    #print "taxid_genomes=", gdt.taxid_genomes
    #print "gdt.taxids=", gdt.taxids
    print >>sys.stderr, "genomes_records length(total seq files)=", len(genomes_records)

    for taxid in gdt.taxids:
      tmp_egl = 0
      for gi in gdt.taxid_genomes[taxid]:
        #print "taxid=", taxid, "gi=", gi
        #tmp_egl += (len(genomes_records[gi].seq)-rdt.read_length) # no read_length correction is preferred for contigs
        tmp_egl += len(genomes_records[gi].seq) # no read_length correction is preferred

      self.egl.append( tmp_egl ) #effective genome length, Hopefully no one is longer than any genome
    """
    
    #for taxid in gdt.taxids:
    #  self.egl.append( gdt.taxid_length[taxid] )   #this info already in gdata file

    ### initiate pMONO
    print >>sys.stderr, "pMONO: self.N=%d, self.M=%d" % ( self.N, self.M )
    self.em = gemcore.cMONO(self.N, self.M, gemcore.VectorFloat(rand_f(self.M)))
    self.em.get_P()
    self.em.get_Prg()
    self.init_Prg()
    self.em.init_ENANMN()

  def em_free(): # related to data structure backend
    self.em.del_P()
    self.em.del_Prg()
    del self.em
    #gemcore.delocate( self.em.P, self.N )
    #gemcore.delocate( self.em.Prg, self.N )

class URAE_EM(pMEM):
  #Unknown Reads Abundance Estimation, in K-tuple
  K = None      # kmer length

  # prepare in working mode, run seed() first
  def prep(self, rdt, gdt, k = 6): 
    """ rdt is an object of Read_Data type, 
        required: rdt.reads_file, rdt.read_length
        gdt is an object of Genome_Data type
        required: gdt.genomes_file, gdt.taxid_genomes
    """

    #subprocess.Popen( "pwd", shell=True )
    #subprocess.Popen( "ls -l", shell=True )
    self.gdt = gdt
    self.K = k
    self.M = len(gdt.taxid_genomes) # Get # of genomes in means of taxid 

    if not os.path.exists( rdt.reads_file ): 
      print rdt.reads_file+"doesn't exist!!!"
    reads_handle = open( rdt.reads_file, 'r' )
    reads_records = []
    for seq_record in SeqIO.parse(reads_handle, "fasta"):
      reads_records.append(seq_record)
      self.rl.append(len(seq_record.seq))  # Get rl from here
    self.N = len(reads_records)            # Get # of reads from reads_file records
    
    print "using genome file:", gdt.genomes_file

    if os.environ.get("TMPDIR") != None:
      print "solving concurrent gdt.genomes_file access problem by using %s in each node" % os.environ.get("TMPDIR") 
      r = subprocess.Popen( "cp " + gdt.genomes_file + " $TMPDIR/genomes.fna", shell=True ).communicate()
      if r != (None, None):
        print "return =", r, ", error!"
        print "Error copying temporary genomes.fna file, quitting..."
        exit()
      gdt.genomes_file = os.environ.get("TMPDIR") + "/genomes.fna"
    print "continuing reading genomes.fna file"

    genomes_handle = open( gdt.genomes_file, "r" )
    genomes_records = []
    for seq_record in SeqIO.parse(genomes_handle, "fasta"):
      genomes_records.append(seq_record)

    crw = gemcore.get_HHM()
    crw.disown() #suppress warnings
      #
    i = 0
    for read in reads_records:
        # VOM lines
        #crw.push_back( gemcore.kmer2omap( str(read.seq), k, False ) )
        # UUM lines
      tm = gemcore.kmer2hmap( str(read.seq), k, False ) 
      tm.disown() # suppress memory leak warnings
        #print tm
        #gemcore.ins_UUM( crw, i, gemcore.kmer2umap( str(read.seq), k, False ) )
      gemcore.ins_HHM( crw, i, tm )
      gemcore.del_HM( tm )
        #
      i += 1
    del reads_records

    pgw = gemcore.allocate(self.M, 4**k)
    i = 0

    #print "taxid_genomes=", gdt.taxid_genomes
    print "gdt.taxids=", gdt.taxids
    print "genomes_records length(total seq files)=", len(genomes_records)

    for taxid in gdt.taxids:
      tmp_gl = 0
      tmp_egl = 0
      tmp_cgw = np.zeros(4**k)
      for gi in gdt.taxid_genomes[taxid]:
        print "taxid=", taxid, "gi=", gi
        tmp_gl += len(genomes_records[gi].seq)
        tmp_egl += 2*(len(genomes_records[gi].seq)-rdt.read_length)
        #print str(genomes_records[gi].seq)
        tmp_cgw += np.array( gemcore.add2vec( \
          gemcore.kmer2vec( str(genomes_records[gi].seq), k, False),\
          gemcore.kmer2vec( str(genomes_records[gi].seq.complement()), k, False) ) ) 
          #pgw.push_back( gemcore.VectorDouble( np.array(tmp_cgw)/tmp_egl ) )
          #print len(tmp_cgw)
      tmp_cgw = [ v/tmp_egl for v in tmp_cgw ]
      tmp_cgw = [ v/(sum(tmp_cgw)+0.1**k) for v in tmp_cgw ]
      for j in range(0, len(tmp_cgw)):
        gemcore.ins_AF( pgw, i, j, tmp_cgw[j] )

        #for j in range(0, len(tmp_cgw)): 
          #print tmp_cgw
          #if tmp_egl == 0: break   #each pgw=0 since this taxid has no genomes associate with it
          #should be some tmp_egl=0 here if there is taxid which has no genomes in ncbi
          #kind of wasting some memory here, but the size of this matrix is not a big issue yet
        #  val = tmp_cgw[j]/tmp_egl
        #  gemcore.ins_AF( pgw, i, j, val )

      self.egl.append( tmp_egl ) #effective genome length, Hopefully no one is longer than any genome
      self.gl.append( tmp_gl )
      i += 1

      #adjust dummy genome size by some constant, here 3.7M as in E.coli
      #if gdt.taxids[-1] == 'dummy':
      #  self.egl[-1] = 7400000
      #handles 'dummy' mode now in solva

      #self.s = v
    
      #old mdmd_URAE way
      #pij = gemcore.MatrixDouble() #consider which creation is better here
      #for i in range(0, self.N):
      #  pij.push_back( gemcore.VectorDouble(self.M, 0.) )
   
      ### final we can initialize sim
      #self.em = gemcore.mdmd_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw) 

      #new adad_URAE way
    pij = gemcore.allocate( self.N, self.M )
    prg = gemcore.allocate( self.N, self.M )
      #VOM lines
      #self.em = gemcore.adadvo_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw)
      #UUM lines
      #print crw
      #print gemcore.adaduu_URAE
      #self.em = gemcore.afafuu_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw)
    self.s = rand_f(self.M) 
    self.em = gemcore.afafhh_URAE(self.K, self.N, self.M, gemcore.VectorFloat(self.s), pij, pgw, crw, prg)
    self.em.init_prg()
    gemcore.del_HHM( crw )
      #
      #gemcore.delocate( pij, self.N )
      #gemcore.delocate( pgw, self.M )
    del pij
    del pgw
      # temporarily here, print gl
    print "taxids=", gdt.taxids

      #get self.gl
      #

    del genomes_records

class URAE_EM(pMEM):
  #Unknown Reads Abundance Estimation, in K-tuple
  K = None      # kmer length

  # prepare in working mode, run seed() first
  def prep(self, rdt, gdt, k = 6): 
    """ rdt is an object of Read_Data type, 
        required: rdt.reads_file, rdt.read_length
        gdt is an object of Genome_Data type
        required: gdt.genomes_file, gdt.taxid_genomes
    """

    #subprocess.Popen( "pwd", shell=True )
    #subprocess.Popen( "ls -l", shell=True )
    self.gdt = gdt
    self.K = k
    self.M = len(gdt.taxid_genomes) # Get # of genomes in means of taxid 
    self.rl = []
    self.egl = []
    self.gl = []

    if not os.path.exists( rdt.reads_file ): 
      print rdt.reads_file+"doesn't exist!!!"
    reads_handle = open( rdt.reads_file, 'r' )
    reads_records = []
    for seq_record in SeqIO.parse(reads_handle, "fasta"):
      reads_records.append(seq_record)
      self.rl.append(len(seq_record.seq))  # Get rl from here
    self.N = len(reads_records)            # Get # of reads from reads_file records
    
    print "using genome file:", gdt.genomes_file

    if os.environ.get("TMPDIR") != None:
      print "solving concurrent gdt.genomes_file access problem by using $TMPDIR=%s in each node" % os.environ.get("TMPDIR") 
      r = subprocess.Popen( "cp " + gdt.genomes_file + " $TMPDIR/genomes.fna", shell=True ).communicate()
      if r != (None, None):
        print "return =", r, ", error!"
        print "Error copying temporary genomes.fna file, quitting..."
        exit()
      gdt.genomes_file = os.environ.get("TMPDIR") + "/genomes.fna"
    print "continuing reading genomes.fna file"

    genomes_handle = open( gdt.genomes_file, "r" )
    genomes_records = []
    for seq_record in SeqIO.parse(genomes_handle, "fasta"):
      genomes_records.append(seq_record)

    crw = gemcore.get_HHM()
    crw.disown() #suppress warnings
      #
    i = 0
    for read in reads_records:
        # VOM lines
        #crw.push_back( gemcore.kmer2omap( str(read.seq), k, False ) )
        # UUM lines
      tm = gemcore.kmer2hmap( str(read.seq), k, False ) 
      tm.disown() # suppress memory leak warnings
        #print tm
        #gemcore.ins_UUM( crw, i, gemcore.kmer2umap( str(read.seq), k, False ) )
      gemcore.ins_HHM( crw, i, tm )
      gemcore.del_HM( tm )
        #
      i += 1
    del reads_records

    pgw = gemcore.allocate(self.M, 4**k)
    i = 0

    #print "taxid_genomes=", gdt.taxid_genomes
    #print "gdt.taxids=", gdt.taxids
    print "genomes_records length(total seq files)=", len(genomes_records)

    for taxid in gdt.taxids:
      tmp_gl = 0
      tmp_egl = 0
      tmp_cgw = np.zeros(4**k)
      for gi in gdt.taxid_genomes[taxid]:
        #print "taxid=", taxid, "gi=", gi
        tmp_gl += len(genomes_records[gi].seq)
        tmp_egl += 2*(len(genomes_records[gi].seq)-rdt.read_length)
        #print str(genomes_records[gi].seq)
        tmp_cgw += np.array( gemcore.add2vec( \
          gemcore.kmer2vec( str(genomes_records[gi].seq), k, False),\
          gemcore.kmer2vec( str(genomes_records[gi].seq.complement()), k, False) ) ) 
          #pgw.push_back( gemcore.VectorDouble( np.array(tmp_cgw)/tmp_egl ) )
          #print len(tmp_cgw)
      tmp_cgw = [ v/tmp_egl for v in tmp_cgw ]
      tmp_cgw = [ v/(sum(tmp_cgw)+0.1**k) for v in tmp_cgw ]
      for j in range(0, len(tmp_cgw)):
        gemcore.ins_AF( pgw, i, j, tmp_cgw[j] )

        #for j in range(0, len(tmp_cgw)): 
          #print tmp_cgw
          #if tmp_egl == 0: break   #each pgw=0 since this taxid has no genomes associate with it
          #should be some tmp_egl=0 here if there is taxid which has no genomes in ncbi
          #kind of wasting some memory here, but the size of this matrix is not a big issue yet
        #  val = tmp_cgw[j]/tmp_egl
        #  gemcore.ins_AF( pgw, i, j, val )

      self.egl.append( tmp_egl ) #effective genome length, Hopefully no one is longer than any genome
      self.gl.append( tmp_gl )
      i += 1

      #adjust dummy genome size by some constant, here 3.7M as in E.coli
      #if gdt.taxids[-1] == 'dummy':
      #  self.egl[-1] = 7400000
      #handles 'dummy' mode now in solva

    self.s = rand_f(self.M) 
      #self.s = v
    
      #old mdmd_URAE way
      #pij = gemcore.MatrixDouble() #consider which creation is better here
      #for i in range(0, self.N):
      #  pij.push_back( gemcore.VectorDouble(self.M, 0.) )
   
      ### final we can initialize sim
      #self.em = gemcore.mdmd_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw) 

      #new adad_URAE way
    pij = gemcore.allocate( self.N, self.M )
    prg = gemcore.allocate( self.N, self.M )
      #VOM lines
      #self.em = gemcore.adadvo_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw)
      #UUM lines
      #print crw
      #print gemcore.adaduu_URAE
      #self.em = gemcore.afafuu_URAE(self.K, self.N, self.M, self.s, pij, pgw, crw)
    self.em = gemcore.afafhh_URAE(self.K, self.N, self.M, gemcore.VectorFloat(self.s), pij, pgw, crw, prg)
    #self.seed()
    self.em.init_prg()
    gemcore.del_HHM( crw )
      #
      #gemcore.delocate( pij, self.N )
      #gemcore.delocate( pgw, self.M )
    del pij
    del pgw
      # temporarily here, print gl
    print "taxids=", gdt.taxids

      #get self.gl
      #

    del genomes_records
    print "self.gl=", self.gl


class DUAL_EM(pMEM):
  K = None
  genomes_records = None
  C_range = None
  read_kvecs = None

  def prep(self, rdt, gdt, k):
    
    print "[STATUS] in prep..."
    self.gdt = gdt
    self.K = k
    self.M = len(self.gdt.taxids)
    self.genomes_records = []
    self.rl = []

    genomes_handle = open( gdt.genomes_file, "r" )
    for seq_record in SeqIO.parse(genomes_handle, "fasta"):
      self.genomes_records.append(seq_record)

    if not os.path.exists( rdt.reads_file ):
      print rdt.reads_file+"doesn't exist!!!"
    reads_handle = open( rdt.reads_file, 'r' )
    self.read_kvecs = gemcore.get_VUM()

    for seq_record in SeqIO.parse(reads_handle, "fasta"):
      #print seq_record.seq.tostring()
      tm = gemcore.kmer2umap( seq_record.seq.tostring(), self.K, True )
      tm.disown()
      gemcore.ins_VUM( self.read_kvecs, tm )
      self.rl.append(len(seq_record.seq.tostring()))  # Get rl from here

    self.N = len(self.rl)            # Get # of reads from reads_file records
    self.rl = np.array(self.rl)

    self.get_gl() #get self.gl and self.egl
    self.s = rand_f(self.M)
    #self.C_range=xrange(self.M,2*self.M+1)
    print "[STATUS] out prep..."
    # after prep, we have all information to find clusters 

  def cluster(self):

    print "[STATUS] in cluster..."

    print "N=", self.N, "M=", self.M
    max_exp = int(np.floor( np.log2( self.N ) ))
    self.C_range=[ 2**e for e in xrange(0, 8) ]
    print "C_range=", self.C_range

    weight = gemcore.VectorFloat( 4**self.K, 1 )
    npass = 10
    dist = 'e'
    error = gemcore.get_flt() #should use extern later
    cids = gemcore.VectorInt( self.N )
    tcids = gemcore.VectorInt( self.N )
    mapping = gemcore.VectorInt( self.N )
    errors = np.zeros( len(self.C_range) )
    aics = np.zeros( len(self.C_range) )
    bics = np.zeros( len(self.C_range) )
    cdatas = []

    print "[STATUS] in clustering..."
    start_time = time.time()
    for i in xrange(0, len(self.C_range)):
      print "round %d" % self.C_range[i]
      #cdata = gemcore.allocate( self.C_range[i], 4**self.K)
      cdata = gemcore.MatrixFloat( self.C_range[i], gemcore.VectorFloat( 4**self.K, 0 ) )  
      counts = gemcore.VectorInt( self.C_range[i] )
      ifound = gemcore.vum_kmeans( self.C_range[i], self.N, 4**self.K, self.read_kvecs, weight, npass, dist, cdata, cids, error, tcids, counts, mapping)
      errors[i] = gemcore.vflt(error)
      aics[i] = gemcore.vflt(error) + 2*self.C_range[i]*(4**self.K)
      bics[i] = self.N * np.log(gemcore.vflt(error)/self.N) + self.C_range[i]*(4**self.K) * np.log(self.N)
      cdatas.append(cdata)
    end_time = time.time()
    print "total time: %d" % (end_time - start_time)
    print "[STATUS] out clustering..."

    print "errors=[", "\t".join([ str(v) for v in errors.tolist() ]), "]"
    print "aics=[", "\t".join([ str(v) for v in aics.tolist() ]), "]"
    print "bics=[", "\t".join([ str(v) for v in bics.tolist() ]), "]"
    Ci = errors.argmin()
    Ai = aics.argmin()
    Bi = bics.argmin()
    C = self.C_range[Ci]
    A = self.C_range[Ai]
    B = self.C_range[Bi]
    print "C=", C, "A=", A, "B=", B
    
    return cdatas, C, Ci

  #make a blackbox here,
  #the input is Cdatas[Ci], C, self.M, self.genomes_records, self.gdt.taxids, self.rl, self.K
  #the output is mcg

  def classify_all(self, cdata, C):  #cdata should be in numpy.array format

    #mcg = gemcore.allocate(C, self.M)
    mcg_ary = np.zeros( (C, self.M) )
    mc = gemcore.VectorFloat( C, 0 )

    for i in xrange(0, self.M): # use 'N' to concatenate the genomes
      print "[STATUS: classifying gi=%d" % i, "len=%d" % self.egl[i] 
      start_time = time.time()
      gemcore.gclassify(C, cdata, \
          'N'.join( [ self.genomes_records[gi].seq.tostring() for gi in self.gdt.taxid_genomes[self.gdt.taxids[i]] ] ), \
          int(math.ceil(np.mean(self.rl))), self.K, mc)
      end_time = time.time()
      print "total time: %d" % (end_time - start_time)
      #gemcore.vec2col( mcg, mc, i )
      mcg_arg[:,c] = mc

    return mcg_arg

  def init(self, cdata, C, mcg ):

    vf = gemcore.VectorFloat( rand_f(self.M) )
    p = gemcore.allocate( self.N, C )
    prg = gemcore.allocate( self.N, self.M )
    self.em = gemcore.afafvu_pMEM( vf, p, prg, self.K, self.M, self.N, C  )
    self.em.V = vf
    self.em.pV = vf
    self.em.R = 10000
    rl = gemcore.VectorFloat( self.rl.tolist() )
    #rl = gemcore.VectorFloat( [ len(s) for s in rs ] )
    #print gemcore.VectorFloat( [ float( self.rl[i] ) for i in xrange(0, len(self.rl)) ] )
    self.em.init_prg( self.read_kvecs, cdata, mcg, rl, self.N, self.M, C, self.K )
    print "[STATUS] out init..."

  def get_gl(self):
    self.egl = np.zeros(self.M)
    self.gl = np.zeros(self.M)
    i = 0
    for taxid in self.gdt.taxids:
      tmp_gl = 0
      tmp_egl = 0
      tmp_cgw = np.zeros(4**self.K)
      for gi in self.gdt.taxid_genomes[taxid]:
        try:
          tmp_gl += len(self.genomes_records[gi].seq.tostring())
        except IndexError:
          print "gi=", gi
          exit()
        tmp_egl += 2*(len(self.genomes_records[gi].seq.tostring())-np.mean(self.rl))
      if taxid == 'dummy':
        tmp_egl = 8000000 # test if this will affect result
        tmp_gl = tmp_egl/2
      self.egl[i] = tmp_egl
      self.gl[i] = tmp_gl
      i += 1

"""
  URAE_EM.prep_sim
  #preparing in simulation mode, run seed() first
  #deprecated code, no use
  def prep_sim(self, gdt, abd, rn, rl, k=6):
    #assume gdt, abd are full, only est is partial, all of them are pre sorted by abd
    self.gdt = gdt
    self.K = k
    self.M = len(self.gdt.taxids)
    genomes_handle = open( gdt.genomes_file, "r" )
    genomes_records = []
    for seq_record in SeqIO.parse(genomes_handle, "fasta"):
      genomes_records.append(seq_record)

    pgw_copy = np.zeros( (self.M, 4**self.K ) )
    pgw = gemcore.MatrixFloat()
    pgw = gemcore.allocate(self.M, 4**self.K)
    i = 0
    for taxid in self.gdt.taxids:
      tmp_gl = 0
      tmp_egl = 0
      tmp_cgw = np.zeros(4**k)
      for gi in gdt.taxid_genomes[taxid]:
        try:
          tmp_gl += len(genomes_records[gi].seq)
        except IndexError:
          print "gi=", gi
          exit()
        tmp_egl += 2*(len(genomes_records[gi].seq)-rl)
        tmp_cgw += np.array( gemcore.add2vec( \
          gemcore.kmer2vec( str(genomes_records[gi].seq), k, False),\
          gemcore.kmer2vec( str(genomes_records[gi].seq.complement()), k, False) ) )
      tmp_cgw = [ v/tmp_egl for v in tmp_cgw ]
      tmp_cgw = [ v/(sum(tmp_cgw)+0.1**k) for v in tmp_cgw ]
      for j in range(0, len(tmp_cgw)):
        gemcore.ins_AF( pgw, i, j, tmp_cgw[j] )
        pgw_copy[i,j] = tmp_cgw[j]
      if taxid == 'dummy':
        tmp_egl = 8000000 # test if this will affect result
        tmp_gl = tmp_egl/2
      self.egl.append( tmp_egl ) 
      self.gl.append( tmp_gl )
      i += 1
    del genomes_records

    rr = [ abd[i]*self.egl[i] for i in xrange(0, self.M) ]
    rr = [ int(r/sum(rr)*rn) for r in rr ]

    self.N = sum(rr)

    crw = gemcore.get_HHM()
    crw.disown()
    cnt = 0

    #print "pgw=", pgw_copy
    #print "crw=",

    for i in xrange(0, self.M):
      #print "pgw_copy[%d]" % i, sum(pgw_copy[i]), pgw_copy[i]
      kvecs = np.random.multinomial( rl, pgw_copy[i], rr[i] ) # prduce rr[i] samples
      for v in kvecs:
        #print "vf=", v
        tm = gemcore.vec2hm( gemcore.VectorFloat(v.tolist()) )
        tm.disown()
        #gemcore.showhmap( tm )
        gemcore.ins_HHM( crw, cnt, tm )
        gemcore.del_HM( tm )
        cnt += 1

    pij = gemcore.allocate( self.N, self.M )
    prg = gemcore.allocate( self.N, self.M )

    self.s = rand_f(self.M) 
    self.em = gemcore.afafhh_URAE(self.K, self.N, self.M, gemcore.VectorFloat(self.s), pij, pgw, crw, prg)
    #run self.seed() or slef.reset() before self.prep_sim()
    self.em.init_prg()
    gemcore.del_HHM( crw )
    #gemcore.delocate( pij, self.N )
    #gemcore.delocate( pgw, self.M )
    del pij
    del pgw
"""
