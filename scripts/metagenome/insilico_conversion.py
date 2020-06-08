# Does insilico conversion ...
import sys
import os
import itertools
from Bio.Seq import Seq

infile=sys.argv[1]
outfile=sys.argv[2]
read = sys.argv[3]

def read_conversion(sequence, read):
	sequence = sequence.upper()
	if read == 'R2':
		endl = ''
		if (sequence[-1]=='\n'):
			endl = '\n'
			sequence = sequence.strip('\n')
		sequence = str(Seq(sequence).reverse_complement()) + endl
	return(sequence)

def fasta(infile, outfile, read):
	with open(infile) as f, open(outfile, 'w') as w:
		for read_id, sequence in itertools.zip_longest(*[f]*2):
			w.write(read_id)
			sequence = read_conversion(sequence, read)
			w.write(sequence.replace('C', 'T'))
def fastq(infile, outfile, read):
	with open(infile) as f, open(outfile, 'w') as w:
		for read_id, sequence, plus, qual in itertools.zip_longest(*[f]*4):
			w.write(read_id)
			sequence = read_conversion(sequence, read)
			w.write(sequence.replace('C', 'T'))
			w.write(plus)
			w.write(qual)
f = open(infile)
line = f.readline()
f.close()
if line[0] == '>':
	fasta(infile, outfile, read)
if line[0] == '@':
	fastq(infile, outfile, read)
