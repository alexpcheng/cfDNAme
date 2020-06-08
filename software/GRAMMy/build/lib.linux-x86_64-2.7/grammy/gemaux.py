#Lsat Updated: Aug 25, 2010
#ADD: init function for Read_Data and Genome_Data type 

import time, os, sys, subprocess, csv, re, random, string
import tarfile, gzip, fnmatch
from Bio import SeqIO

#global rdata, gdata #reserved for rdt and gdt file import

class Read_Data:
  """base read data class"""
  def __init__(self):
    self.read_length = None      # average read length
    self.read_tech = None        # sequencing technology
    self.reads_file = None       # file contains reads
    self.reads_number = None     # either use the # of seqs from reads_file or this
    self.reads_files = None      # if the reads data is a set of reads file
  
  def write(self, prefix):
    rdt_file = open(prefix+".rdt", 'w')
    print >>rdt_file, "rdata=Read_Data()"
    print >>rdt_file, "rdata.read_tech= '"+str(self.read_tech)+"'"
    print >>rdt_file, "rdata.reads_file= '"+str(self.reads_file)+"'"
    print >>rdt_file, "rdata.read_length= %d" % self.read_length
    print >>rdt_file, "rdata.reads_number= %d" % self.reads_number

  def read(self, prefix):
    exec open(prefix+".rdt").read()
    self.read_length = rdata.read_length
    self.read_tech = rdata.read_tech
    self.reads_file = rdata.reads_file
    self.reads_number = rdata.reads_number

class Genome_Data:
  """base genome data class"""
  def __init__(self):
    self.taxids = None           # taxids of interest
    self.taxid_genomes = None    # taxid_genomes correspondence 
                                 # one taxid may corespondence to several genomes and each genome with lots of sequence 
    self.gi_taxid = None         # gi_taxid correspondence, each gi has one taxid but each taxid may have a bunch of gis
    self.genomes_file = []       # file containing sequences of all genomes, can be in series from .1, .2, ...
    self.taxids_number = None    # either use this or use len(taxids) or len(taxid_genomes)
    self.taxid_length = None     # genome length for taxid

  def write(self, prefix):
    gdt_file = open(prefix+".gdt", 'w')
    print >>gdt_file, "gdata=Genome_Data()"
    print >>gdt_file, "gdata.taxids=" +str(self.taxids)
    print >>gdt_file, "gdata.taxid_genomes= "+str(self.taxid_genomes)
    print >>gdt_file, "gdata.gi_taxid=" +str(self.gi_taxid)
    print >>gdt_file, "gdata.genomes_file="+str(self.genomes_file)
    print >>gdt_file, "gdata.taxid_length="+str(self.taxid_length)

  def read(self, prefix):
    exec open(prefix+".gdt").read()
    self.taxids = gdata.taxids
    self.taxid_genomes = gdata.taxid_genomes
    self.gi_taxid = gdata.gi_taxid
    self.genomes_file = gdata.genomes_file
    self.taxids_number = gdata.taxids_number
    self.taxid_length = gdata.taxid_length

class Sim_Read_Data(Read_Data):  # inherit not for memebers
  """sim read data class"""
  def __init__(self):
    #self.read_length = None      # average read length
    #self.read_tech = None        # sequencing technology
    #self.reads_file = []         # file contains reads
    #self.reads_number = None     # either use the # of seqs from reads_file or this
    #self.reads_files = None      # if the reads data is a set of reads file
    self.gdt = Genome_Data()     # Genome_Data used to simulate the reads
                                 # require: gdt.taxids, gdt.taxid_genomes, gdt.gi_taxid, gdt.genomes_file, gdt.genome_number
    self.taxid_reads = None      # taxid_reads correspondence of simulated reads
    self.abundance = None        # abundance vector used to simulate the reads
    self.reads_numbers = None    # reads_number for each taxid, different from reads_number
    self.error_model = None      # error model used to generate the reads set

def Read_by_File( rf, rt, rl, format='fasta' ):
  rdt = Read_Data()
  rdt.read_tech = rt
  rdt.read_length = rl
  rdt.reads_file = rf
  rn = 1
  for seq in SeqIO.parse( open(rf, 'r'), format):
    rn += 1
  rdt.reads_number = rn
  return rdt

def Genome_by_gref( taxids, o_prefix, dmp_file="grefs/gid_tid.dmp", gref_dir="grefs", per_set=20 ):
  gdt = Genome_Data()
  gi_2_ti = dict()
  taxids_set = set(taxids)

  for row in open( dmp_file ):
    gid, tid = row.rstrip('\n').split('\t')
    if int(tid) in taxids_set:
      gi_2_ti[int(gid)] = int(tid)

  genomes_out_tpl = o_prefix + ".fna.%d"
  genomes_file = []
  ti_gl = dict()
  ti_2_record = dict()
  cnt = 0
  for i in range(0, len(taxids)):
    if i % per_set == 0:
      #close output file and renew another
      genomes_out_handle = open( genomes_out_tpl % ( i/per_set+1 ), 'w' )
      genomes_file.append( genomes_out_tpl % ( i/per_set+1 ) )
    assert os.path.exists( os.path.join( gref_dir, str(taxids[i]), 'genomes.fna' ) )
    for seq_rec in SeqIO.parse( open( os.path.join( gref_dir, str(taxids[i]), 'genomes.fna' ) ), "fasta" ):  #taxid must be in grefs
      ti_2_record.setdefault( taxids[i], [] ).append( cnt )
      ti_gl[taxids[i]] = ti_gl.get(taxids[i],0) + len(seq_rec)
      cnt += 1
      SeqIO.write([seq_rec], genomes_out_handle, 'fasta')
  print >>sys.stderr, "loaded %d" % cnt, "seqs into %d" % (len(taxids)/per_set+1), "files"
      
  gdt.taxid_genomes = ti_2_record
  gdt.gi_taxid = gi_2_ti
  gdt.taxids = taxids
  gdt.taxid_length = ti_gl
  gdt.taxid_number = len(taxids)
  gdt.genomes_file = genomes_file
  return gdt

def Genome_by_NCBI( taxids, o_prefix, taxids_universe, dummy_mode=False, dmp_file="gi_taxid_nucl.dmp", fna_file=None, add_file=None,
    complete_dir='genbank/Bacteria', draft_dir=None ): 

  # random string to avoid parallel collision
  tmpdir = os.environ.get('TMPDIR')
  print >>sys.stderr, "---- tmpdir %s" % tmpdir

  # constructing Genome_Data from standard files of NCBI taxonomy
  # drop those ids that do not 
  gdt = Genome_Data()

  # genomes output file name
  genomes_out_name = o_prefix+".fna"

  # read gi 2 ti dump and create ti_2_gi map
  ti_2_gen=dict()
  for ti in taxids:
    ti_2_gen[ti] = []
  gi_2_ti = dict()

  taxids_set = set( taxids ) # for fast taxid index
  #print "reading in dmp file, this may take a while"
  print >>sys.stderr, "Getting for taxid <-> genome correspondence..."
  dmp_handle = None
  start_time = time.time()
  if dmp_file.endswith(".gz"):
    args = "gzip -dc %s" % dmp_file
    dmp_handle = subprocess.Popen(args, stdout=subprocess.PIPE, shell=True).stdout
  else:
    dmp_handle = open( dmp_file, "r" )
  for row in dmp_handle:
    gid, tid = row.rstrip('\n').split('\t')
    if int(tid) in taxids_set:
      gi_2_ti[int(gid)] = int(tid) #save memory
    else:
      pass
  dmp_handle.close()
  end_time = time.time()
  print >>sys.stderr, "done finding correspondence"
  print >>sys.stderr, "Time: \n %d s" % ( end_time - start_time )
  
  # read genome database and create gi_2_record map
  print >>sys.stderr, "Looking for genomes..."
  start_time = time.time()
  genomes_out_handle = open( genomes_out_name, "w" )
  #genomes_records = []
  ti_2_record = dict()
  ti_gl = dict()

  # STEP(1): add gzipped genomes
  rkey = ''.join( random.sample( string.letters, 8 ) )
  print >>sys.stderr, "---- rkey_one %s" % rkey
  tid = None
  gid = None
  ct = 0
  if fna_file != None:
    args = 'rm -rf '+tmpdir+'/'+rkey
    tar = tarfile.open( fna_file )
    tar.extractall( path=tmpdir+'/'+rkey )
    tar.close()
    #avoid big genomes_records data here
    for root, dirs, files in os.walk(tmpdir+'/'+rkey):
      for name in files:
        if name.endswith('.fna') or name.endswith('.fsa_nt'):  # .fsa_nt for contigs
          #print >>sys.stderr, "\tParsing %s" % name
          tmp_fnafile = open( os.path.join( root, name ), "r" )
          for seq_rec in SeqIO.parse( tmp_fnafile , 'fasta'):
            gid = int( seq_rec.description.split('|')[1] )
          if gid in gi_2_ti:
            tid = gi_2_ti[gid]
            ti_2_record.setdefault( tid, [] ).append(ct)
            ti_gl[tid] = ti_gl.setdefault( tid, 0 ) + len(seq_rec) 
            ct += 1
            SeqIO.write([seq_rec], genomes_out_handle, 'fasta')
          tmp_fnafile.close()
  print >>sys.stderr, "loaded %d" % ct, "seqs from compressed fna file"

  # STEP(1), (2), (3) sequence should be exchanged?

  # STEP(2): add complete genomes 
  tid = None
  gid = None
  #print ti_2_record
  rkey = ''.join( random.sample( string.letters, 8 ) )
  print >>sys.stderr, "---- rkey_two = %s" % rkey
  print >>sys.stderr, "---- complete_dir = %s" % complete_dir
  bt = ct
  if complete_dir != None:
    for subdir in os.listdir( complete_dir ):
      tmp_seqs = []
      tid = None
      tmp_tid = None
      tmp_gid = None
      bail = False
      for file in os.listdir( os.path.join( complete_dir, subdir ) ):
        if file.endswith('.fna'):
          tmp_fnafile = open( os.path.join( complete_dir, subdir, file ) ) 
          for seq_rec in SeqIO.parse( tmp_fnafile, 'fasta' ):
            tmp_gid = int(seq_rec.description.split('|')[1])
            try:
              tmp_tid = gi_2_ti[tmp_gid]
            except KeyError:
              bail = True # quit to next dir
            assert tid == None or tid == tmp_tid  #ensure all seqs come from the same tid
            tid = tmp_tid
            tmp_seqs.append( seq_rec )
          tmp_fnafile.close()
      if tid in ti_2_record.keys():
        print >>sys.stderr, "skipping duplicated record %s" % ( tid )
        continue
      elif tid == None:
        print >>sys.stderr, "skipping missing gids \n%s" % ( '\n'.join( [ seq_rec.description for seq_rec in tmp_seqs ] ) )
        continue
      elif bail:
        print >>sys.stderr, "partial missing gids \n%s" % ( '\n'.join( [ seq_rec.description for seq_rec in tmp_seqs ] ) )
        continue
      else:
        for seq_rec in tmp_seqs:
          #tmp_gid = int(seq_rec.description.split('|')[1])
          ti_2_record.setdefault( tid, [] ).append(bt)
          ti_gl[tid] = ti_gl.setdefault( tid, 0 ) + len(seq_rec)
          bt += 1
          SeqIO.write([seq_rec], genomes_out_handle, 'fasta')
          #print >>sys.stderr, 'added: gi=%d, ti=%d' % (tmp_gid, tid)
      print >>sys.stderr, "subdir %s, updated:" % subdir
  print >>sys.stderr, "loaded %d" % bt, "seqs from genbank/Bacteria"

  #STEP(3): add draft genomes
  tid = None
  gid = None
  rkey = ''.join( random.sample( string.letters, 8 ) )
  print >>sys.stderr, "---- rkey_three = %s" % rkey
  print >>sys.stderr, "---- draft_dir = %s" % draft_dir
  dt = bt
  if draft_dir != None:
    for subdir in os.listdir( draft_dir ):
      tid, tmp_tid, tmp_gid, bail, tmp_seqs = None, None, None, False, []
      for file in os.listdir( os.path.join( draft_dir, subdir ) ):
        if file.endswith('.contig.fna.tgz'):
          if os.path.exists( os.path.join( tmpdir, rkey ) ):
            args = 'rm -rf '+ os.path.join( tmpdir, rkey )
            print >>sys.stderr, "Sys Call:" + args
            subprocess.Popen( args, shell=True ).communicate()
          tmp_tar = tarfile.open( os.path.join( draft_dir, subdir, file ) ) 
          tmp_tar.extractall( path=os.path.join( tmpdir, rkey ) )
          for root, dirs, files in os.walk( os.path.join( tmpdir, rkey ) ):
            for file in files:
              tmp_fnafile = open( os.path.join( tmpdir, rkey, file ) ) 
              for seq_rec in SeqIO.parse( tmp_fnafile, 'fasta' ):
                tmp_gid = int( seq_rec.description.split('|')[1] )
                try:
                  tmp_tid = gi_2_ti[tmp_gid]
                except KeyError:
                  bail = True
                assert tid == None or tid == tmp_tid  #ensure all seqs come from the same tid
                tid = tmp_tid
                tmp_seqs.append( seq_rec )
              tmp_fnafile.close()
      if tid in ti_2_record.keys():
        print >>sys.stderr, "skipping duplicated record %s" % ( tid )
        continue
      elif tid == None:
        print >>sys.stdout, "skipping missing gids \n%s" % ( '\n'.join( [ seq_rec.description for seq_rec in tmp_seqs ] ) )
        continue
      elif bail:
        print >>sys.stdout, "skipping partial missing gids \n%s" % ( '\n'.join( [ seq_rec.description for seq_rec in tmp_seqs ] ) )
        continue
      else:
        for seq_rec in tmp_seqs:
          #tmp_gid = seq_rec.description.split('|')[1]
          ti_2_record.setdefault( tid, [] ).append(dt)
          ti_gl[tid] = ti_gl.setdefault( tid, 0 ) + len(seq_rec)
          dt += 1
          SeqIO.write([seq_rec], genomes_out_handle, 'fasta')
          #print >>sys.stderr, 'added: gi=%s, ti=%s' % (tmp_gid, tid)
      print >>sys.stderr, "subdir %s, updated:" % subdir
  print >>sys.stderr, "loaded %d" % (dt-bt), "seqs from genbank/Bacteria_DRAFT"

  #STEP(4): add customized genomes
  #following need to be modified
  tid = None
  gid = None
  if add_file != None: # add contigs or self prepared genomes, ensure one taxid is in one sub_dir
    rkey = ''.join( random.sample( string.letters, 8 ) )
    print >>sys.stderr, "---- rkey_four = %s" % rkey
    print >>sys.stderr, "---- add_file = %s" % add_file
    et = dt
    print >>sys.stderr, "adding customized genomes..."
    add_tar = tarfile.open( add_file )
    add_tar.extractall( path=os.path.join( tmpdir,rkey ) )
    add_tar.close()
    # add_file format taxid/seqs
    for subdir in os.listdir( os.path.join( tmpdir, rkey ) ):
      tmp_seqs = []
      tid = int( subdir )
      for file in os.listdir( os.path.join( tmpdir, rkey, subdir ) ): 
        if file.endswith('.fna') or name.endswith('.fsa_nt'):
          tmp_fnafile = open( os.path.join( tmpdir, rkey, subdir, file ) )
          for seq_rec in SeqIO.parse( tmp_fnafile ):
            tmp_seqs.append( seq_rec )
          tmp_fnafile.close()
      if tid in ti_2_record.keys():
        print >>sys.stderr, "skipping duplicated record %s" % ( tid )
        continue
      else:
        for seq_rec in tmp_seqs:
          #tmp_gid = seq_rec.description.split('|')[1]
          ti_2_record.setdefault( tid, [] ).append(et)
          ti_gl[tid] = ti_gl.setdefault( tid, 0 ) + len(seq_rec)
          et += 1
          SeqIO.write([seq_rec], genomes_out_handle, 'fasta')
    print >>sys.stderr, "loaded %d" % (et-dt), "seqs from add_file"

  genomes_out_handle.close()
  end_time = time.time()
  #print taxids
  #print ti_2_record
  #exit(-1)
  print >>sys.stderr, "Time: \n %d s" % ( end_time - start_time )
  #print "=========EOSC========="

  gdt.genomes_file = genomes_out_name  
  gdt.taxid_genomes = ti_2_record
  gdt.gi_taxid = gi_2_ti
  #gdt.taxids = gdt.taxid_genomes.keys() # order of taxid_genomes doesn't matter?
  gdt.taxids = taxids
  gdt.taxids_number = len(gdt.taxids)
  gdt.taxid_length = ti_gl
  #gdt.taxids = taxids
  #gdt.taxids_number = len(taxids)
  #if len(gdt.taxids) != len(gdt.taxid_genomes):
  #  print >>sys.stderr, gdt.taxids
  #  print >>sys.stderr, gdt.taxid_genomes.keys()
  #  print >>sys.stderr, "Some ti not found in NCBI, leave them out..."
  return gdt

def Genome_Data_Clean( gdt ):
  """ clean gdt and remove taxids that don't have genomes represented in NCBI refseq """
  present_ids = []
  for ti in gdt.taxids:
    if ti in gdt.taxid_genomes.keys():
      present_ids.append(ti)
  gdt.taxids = present_ids
  gdt.taxids_number = len(gdt.taxids)
  return gdt

# help function to run the simulator: MetaSim
def run_metasim( args ):
  print "Sys Call:" + args
  proc = subprocess.Popen(args, shell=True)
  proc.wait()
  print "metasim ok"

def Read_by_MetaSim( gdt, abd, o_prefix, o_set, rnumber=100, rlength=400, rtech='454', metasim_cmd="$HOME/app/metasim/MetaSim", error=False ):          
  # suppose the MetaSim has already been updated

  out_dir = os.path.split(o_prefix)[0]

  rdt = Sim_Read_Data() 
  rdt.gdt = gdt
  #rdt.taxids = gdt.taxids
  #rdt.genomes_file = gdt.genomes_file
  #rdt.taxid_genomes = gdt.taxid_genomes
  rdt.taxid_reads = dict()
  rdt.read_length = rlength
  rdt.reads_number = rnumber
  rdt.read_tech = rtech
  rdt.abundance = abd 
  rdt.error_model = 'Exact' 
  out_flag = '-Exact'
  error_flag = ''
  if error:
    if rdt.read_length >= 400:
      rdt.error_model = 'Sanger'
      out_flag = '-Sanger'
      error_flag = ' --sanger '
    elif rdt.read_length >= 50:
      rdt.error_model = '454'
      out_flag = '-454'
      error_flag = ' --454 '
    else:
      rdt.error_model = 'Empirical'
      out_flag = '-Empirical'
      error_flag = ' --empirical '
  mprf_out_name = o_prefix+"_L"+str(rlength)+"_N"+str(rnumber)+"_S"+str(o_set)+".mprf"
  reads_in_name = o_prefix+"_L"+str(rlength)+"_N"+str(rnumber)+"_S"+str(o_set)+".fasta"
  metasim_trace = o_prefix+"_L"+str(rlength)+"_N"+str(rnumber)+"_S"+str(o_set)+"%s.tarce.gz" % out_flag
  metasim_out_name = o_prefix+"_L"+str(rlength)+"_N"+str(rnumber)+"_S"+str(o_set)+"%s.fna" % out_flag
  rdt.reads_file = reads_in_name

  mprf_out_file = open( mprf_out_name, "w" )
  for i in range(0,len(gdt.taxids)):
    print >>mprf_out_file, rdt.abundance[i], '\t', "taxid", '\t', gdt.taxids[i]
  mprf_out_file.close() 
  # this must be the reason cause metasim to fail
  # mprf_out_file cannot read by java if it is still open with python
  # a good programming style is open file just before it is needed and close it right after use

  args0 = " rm -f %s %s %s " % (metasim_out_name, reads_in_name, metasim_trace)
  proc0 = subprocess.Popen(args0, shell=True)
  proc0.wait()
  args1 = metasim_cmd + " cmd -s --threads 0 -d %s -r%s -f%s -t0 %s %s 2>&1 >%s " % ( 
      out_dir, str(rdt.reads_number), str(rdt.read_length), error_flag, mprf_out_name, reads_in_name+".log" )
  print >>sys.stderr, args1
  run_metasim( args1 ) 
  args2 = " mv %s %s " % (metasim_out_name, reads_in_name)
  proc2 = subprocess.Popen(args2, shell=True)
  proc2.wait()
  print "run_metasim shall finish before this"

  time.sleep(5)

  args3 = "  rm -f %s " % (metasim_trace)
  #print args3
  proc3 = subprocess.Popen(args3, shell=True)
  proc3.wait()
  
  reads_in_file = open( reads_in_name, "r" )
  n = 0
  for seq_rec in SeqIO.parse(reads_in_file, "fasta"):
    gi = int((re.search('GI=(\d+)', seq_rec.description)).group(1))
    if gdt.gi_taxid[ gi ] in rdt.taxid_reads:
      rdt.taxid_reads[ gdt.gi_taxid[ gi ] ].append(n)
    else:
      rdt.taxid_reads[ gdt.gi_taxid[ gi ] ] = [n] 
    n += 1
  reads_in_file.close()

  return rdt
  # generate genome file taxids only

def Taxid_by_GOLD( keyword='Bacteria|Archaea', ci=0, gold_file="gold_published_table.txt" ):
  """ Pick out the taxids which has keyword in column[ci] in the gold_table.
      the keyworld can be a regular expression such as 'Bacteria|Archaea'.
      It will by default pick out all the Bacteria and Archaea taxids.
      Intersect or Union the returned result to get more control on the taxids you select.
      Index for some important columns in GOLD table:
      'SUPERKINDOM': row[0],
      'PHYLUM' : row[1],
      'CLASS' : row[2],
      'ORDER' : row[3],
      'FAMILY' : row[4],
      'GENUS' : row[5],
      'SPECIES': row[6],
      'Organism': row[8],
      'STRAIN': row[9],
      'SEROVAR': row[10],
      'RELEVANCE': row[21],
      'DISEASE': row[22],
      'HABITAT': row[23],
      'ISOLATION': row[34],
      'TAXID': row[42],
  """
  gold_table = open( gold_file, "r" )
  table_reader = csv.reader(gold_table, delimiter='\t')
  Gold_TaxonId = []
  for row in table_reader:
    if re.search( keyword, row[ci] ):
      Gold_TaxonId.append(row[42])
  return Gold_TaxonId

def Phylo_by_GOLD( taxids, abd, o_prefix, gold_file="gold_published_table.txt" ):
  #Output Phylogenetic tree to a Newick format file
  
  gold_table = open( gold_file, "r" )
  table_reader = csv.reader(gold_table, delimiter='\t')
  Gold_TaxonId = []
  Gold_Info = dict()
  for row in table_reader:
    if not re.search("Error!!!", row[0]):
      #print row[43]
      Gold_TaxonId.append(row[42])
      Gold_Info[ row[42] ] = {  'SUPERKINDOM': row[0],
                                'PHYLUM' : row[1],
                                'CLASS' : row[2],
                                'ORDER' : row[3],
                                'FAMILY' : row[4],
                                'GENUS' : row[5],
                                'SPECIES': row[6],
                                'Organism': row[8],
                                'STRAIN': row[9],
                                'SEROVAR': row[10]
                              } 
  gold_table.close()
  #Format the tree
  infos = [ 'Organism', 'STRAIN', 'SEROVAR' ]
  levels = [ 'SUPERKINDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES', 'STRAIN' ]
  #doc = Document()
  phylo_tree = dict()
  for i in range(0, len(taxids)):
    sub = phylo_tree
    for l in levels:
      if not sub.has_key(Gold_Info[ taxids[i] ][ l ]):
        sub[ Gold_Info[ taxids[i] ][ l ] ] = dict()
        sub = sub[ Gold_Info[ taxids[i] ][ l ] ]
      else:
        sub = sub[ Gold_Info[ taxids[i] ][ l ] ]
      if l == 'STRAIN': 
        sub['Organism'] = ' '.join( [ Gold_Info[ taxids[i] ][ k ] for k in infos ] + \
          [str( abd[i] )] )
  newick = str(phylo_tree)

  def dashrepl( matchobj ):
    if matchobj.group(0) == '{': return '('
    else: return ')'

  newick = re.sub( '[\{|\}]', dashrepl, newick)
  newick = re.sub( '\'[^\']+\':', '', newick)

  NEWICK_FILE = open( o_prefix+".nwk", "w")
  print >>NEWICK_FILE, newick+";" # need to be viewed by Dendroscope
  NEWICK_FILE.close()

  TAXON_FILE = open( o_prefix+".tax", "w")

  #print  "\t".join( levels.append('Organism') ) 
  print >>TAXON_FILE, "\t".join( levels + ['Organism'] )
  for i in range(0, len(taxids)):
    print >>TAXON_FILE, "\t".join( [Gold_Info[ taxids[i] ][ l ] for l in levels] + \
                              [' '.join([ Gold_Info[ taxids[i] ][ k ] for k in infos ]) ] )

  TAXON_FILE.close()


def Fasta_by_Trace( fasta_file, anc_file, trace_file, kw ): 
  """ Reads data from trace fasta and anc file
  """

  def parse_anc( anc_reader ):
    anc_table = []
    #keys = anc_reader.next() 
    # using dictionary passing method will cause MemoryError in python  
    # fall back to list passing
    anc_reader.next()
    #print keys
    for row in anc_reader:
      anc_table.append(row)
      #print anc_table[i]
      #if i>=1: break
    return anc_table

  def pick_seqs( anc_table, trace_reads, kw ):
    picked = []
    i = 0
    for item in anc_table:
      for key in kw.keys():
        if item[key] == kw[key]:
          picked.append(i)
      i += 1
    return [ anc_table[i] for i in picked  ], [ trace_reads[i] for i in picked ] 

  def trim_seqs( anc_table, trace_reads ):
    for i in range(0, len(trace_reads)):
      trim_left = 0
      right_right = len(trace_reads[i].seq) 
      # for details on how to trim see RFC of Trace DB
      # column x=14 is 'CLIP_VECTOR_LEFT'
      # column y=15 is 'CLIP_VECTOR_RIGHT'
      if anc_table[i][14] != '':
        trim_left = int(anc_table[i][14]) - 1 # convert to C notation
      if anc_table[i][15] != '':
        trim_right = int(anc_table[i][15])
      trace_reads[i].seq = trace_reads[i].seq[trim_left:trim_right] # in place replacement 
    return trace_reads
  
  trace_reader = SeqIO.parse( open( trace_file, "r" ), format='fasta' )  
  trace_reads = []
  for read in trace_reader:
    trace_reads.append(read)

  anc_reader = csv.reader( open( anc_file, "r" ), delimiter='\t' )
  anc_table = parse_anc( anc_reader )
  #print anc_table
  #exit()
  if len( anc_table ) != len( trace_reads ):
    print >>sys.stderr, "error!!! not correct pair of anc and fasta file"

  anc_table, trace_reads = pick_seqs( anc_table, trace_reads, kw )

  trace_reads = trim_seqs( anc_table, trace_reads )
  
  fasta_handle = open( fasta_file, "w" )
  SeqIO.write( trace_reads, fasta_handle, "fasta" )
  fasta_handle.close()

class Taxid_Tree:
  #rank_codes = dict({'no rank':0, 'superkindom':1, 'kingdom':2, 'pylum':3, 'class':4, 'order':5, 'family':6, 'ge 
  tree = None

  def __init__( self ):
    self.tree = dict()

  def Tree_from_Nodes( self, dmp_file ):
    pnode_edges = dict()
    for line in dmp_file:
      taxid, parent = line.rstrip('\n').split('|')[0:2]
      taxid = taxid.strip()
      parent = parent.strip()
      print "|"+parent+"|"+"->"+"|"+taxid+"|"
      if parent == "": # root 
        self.tree = { taxid:[] }
      else:
        pnode_edges.setdefault(parent, []).append(taxid)
    print "pnode_edges=",pnode_edges
    #nodes.sort().reverse() # start from leaves nodes and bottom up
  
    def Add_Nodes( sub, pnode_edges ):
      #if not isinstance( sub, dict ): return #meet a leave
      pnode_key = sub.keys()[0]
      print "pnode_key=",pnode_key
      if pnode_key in pnode_edges:
        for cnode_key in pnode_edges[ pnode_key ]:
          sub.setdefault(pnode_key,[]).append({cnode_key:[]})
          #print sub[pnode_key]
          i = len(sub[pnode_key]) - 1
          #print "i=",i
          Add_Nodes( sub[pnode_key][i], pnode_edges )
          #Add_Nodes( sub[cnode_key], pnode_edges )
      else:
        sub[pnode_key]=[0.] # set leaves

    Add_Nodes( self.tree, pnode_edges )
    return self.tree

  def Adjust_Leave_Weight( self, taxid, val ):
    """ Adjust a Leave's Weight through a BFS """
    #pirint "self.tree", self.tree
    stack = [ self.tree ] 
    flag = None
    while len(stack) != 0:
      node = stack.pop()
      #print "node:", node
      key = node.keys()[0]
      #print "key:", key
      if isinstance( node[key][0],float ) or isinstance( node[key][0],int ):
        #  found a leave 
        if key == taxid: # found the correct leave
          node[key][0] = val
          flag = node
          break
      else:
        for child in node[key]:
          #print "child:", child
          stack += [ child ]
          #print "stack:", stack
    return flag

  def Adjust_Node_Weight( self, taxid, val ):
    """ Ajust the weight but don't care the node is leaf or not """
    stack = [ self.tree ] 
    flag = None
    while len(stack) != 0:
      node = stack.pop()
      #print node
      key = node.keys()[0]
      if isinstance( node[key][0],float ) or isinstance( node[key][0],int ):
        if key == taxid:
          node[key] = [val]
          flag = node
          break
      else:
        if key == taxid:
          node[key] = [val]
          flag = node
          break
        for child in node[key]:
          stack += [ child ]
    return flag

  def Find_Branch( self, taxid ):
    """ Find and return a Node through a BFS """
    #pirint "self.tree", self.tree
    stack = [ self.tree ] 
    flag = None
    while len(stack) != 0:
      node = stack.pop()
      #print "node:", node
      key = node.keys()[0]
      #print "key:", key
      if key == taxid: # found the correct leave
        return node
      elif isinstance( node[key][0],float ) or isinstance( node[key][0],int ):
        continue
      else:
        for child in node[key]:
          stack += [ child ]

  def Get_Branch_Weight( self, sub ):
    val = 0
    key = sub.keys()[0]
    if isinstance( sub[key][0],float ) or isinstance( sub[key][0], int ):
      val += sub[key][0]
    else:
      for child in sub[key]:
        val += self.Get_Branch_Weight( child )
    return val

  def Get_Depth_Branches( self, sub, depth, bs ):
    """ get all branches at certain depth by a list """
    print sub
    key = sub.keys()[0]
    if depth == 0:
      #print "depth 0, before bs=:", bs
      bs += [sub]
      #print "depth 0, after bs=:", bs
    else:
      depth -= 1
      #print "sub.values=", sub.values()
      if isinstance( sub[key][0],float ) or isinstance( sub[key][0], int ):
        return #nothing to get any more
      else:
        for child in sub[key]:
          #  print "child=", child
          self.Get_Depth_Branches( child, depth, bs )
          #print "depth=%d" % (depth), branches

  def Get_Depth_Weight( self, sub, depth ):
    branches = []
    self.Get_Depth_Branches( sub, depth, branches )
    vals = []
    for branch in branches:
      vals.append( self.Get_Branch_Weight( branch ) )
    taxids = [ b.keys()[0] for b in branches ]
    return taxids, vals

  def Get_Depth_Tree( self, sub, depth ):
    tree = Taxid_Tree()
    tree.tree = sub
    packed = tree.Get_Depth_Weight( tree.tree, depth ) 
    for i in range(0,len(packed)):
      print packed[0][i], packed[1][i]
      tree.Adjust_Node_Weight( packed[0][i], packed[1][i] )
    return tree

  def Tree_to_NWK( self, sub, depth ):
    tree = self.Get_Depth_Tree( sub, depth )
    newick = str( tree.tree )

    print newick
    
    def dashrepl( matchobj ):
      if matchobj.group(0) == '{': return '('
      else: return ')'

    newick = re.sub( '[\{|\}]', dashrepl, newick)
    #newick = newick.reverse()
    print newick
    newick = re.sub( '\'[^\']+\':', '', newick)

    print newick
