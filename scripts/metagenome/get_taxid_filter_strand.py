import sys
import subprocess
outfile = sys.argv[1]
gi_to_tax = sys.argv[2]
conversion = sys.argv[3]
def get_gi(sseqid):
    gi = sseqid.split('|')[1]
    return(gi)

def get_taxid(gi, gi_to_tax):
    if gi not in gi_to_tax:
        #print(gi)
        taxid = 'GI_NOT_FOUND'
    else:
        taxid = gi_to_tax[gi]
    return(taxid)

with open(gi_to_tax) as f:
    gi_to_tax = dict(x.rstrip().split('\t') for x in f)

with open(outfile, 'w') as w:
    for line in sys.stdin:
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart,qend, sstart, send, evalue, bitscore, qlen, strand = line.strip().split('\t')
        gi = get_gi(sseqid)
        taxid = get_taxid(gi, gi_to_tax)
        if conversion == 'CT' and strand == 'Plus/Plus':
                w.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
        if conversion == 'GA' and strand == 'Plus/Minus':
                w.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
