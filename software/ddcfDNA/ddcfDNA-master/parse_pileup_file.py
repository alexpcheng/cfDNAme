#!/usr/bin/env python3
# Eilon Sharon Oct 2014


"""
###############################################################################
# parses pileup base string and returns the counts for all possible alleles 
# for each position
# reads input (mpileup output) from sys.stdin
# modified from https://github.com/noushin6/NGSTools/blob/master/baseParser.py
# such that it doesn't count the allele when there is an insertion or deletion
###############################################################################
"""

#from __future__ import division



import os
import sys

class parseString(object):
  
  def __init__(self, ref, string):
    self.ref = ref.upper()
    self.string = string.upper()
    self.types = {'A':0,'G':0,'C':0,'T':0,'NotSNP':0,'SNP_Reads':0,'Reads':0,'-':[],'*':0,'+':[],'X':[]}
    self.process()
      
  def process(self):
    # remove end of read character
    self.string = self.string.replace('$','')
    
    
    while self.string != '':
      #sys.stderr.write("DEBUG parse string: %s\n" % self.string)
      
      # 
      if self.string[0] == '^':
        # skip two characters when encountering '^' as it indicates
        # a read start mark and the read mapping quality
        self.string = self.string[2:]
        
      elif self.string[0] == '*':
        self.types['*'] += 1
        # skip to next character
        self.string = self.string[1:]
        self.types['Reads'] += 1
        self.types['NotSNP'] += 1
    
      elif self.string[0] in ['.',',','A','G','C','T']:
        
        self.types['Reads'] += 1

        if (len(self.string)== 1) or (self.string[1] not in ['+','-']):
            if self.string[0] in ['.',',']:
              # a reference base
              self.types[self.ref] += 1
              self.string = self.string[1:]
            else:
              # one of the four bases
              self.types[self.string[0]] += 1
              self.string = self.string[1:]
            
            self.types['SNP_Reads'] += 1

        elif self.string[1] == '+': 
            insertionLength = int(self.string[2])
            insertionSeq = self.string[0] + '+' + self.string[3:3+ insertionLength]
            self.types['+'].append(insertionSeq)
            self.string = self.string[3+insertionLength:]
            self.types['NotSNP'] += 1
            # add cnt for matching nt in insertion if you wish

        elif self.string[1] == '-':
            deletionLength = int(self.string[2])
            deletionSeq = self.string[0] + '-' +  self.string[3:3+deletionLength]
            self.types['-'].append(deletionSeq)
            self.string = self.string[3+deletionLength:]
            self.types['NotSNP'] += 1
            # add cnt for matching nt in deletion if you wish
            
      #elif self.types.has_key(self.string[0]) and\
      #     ((len(self.string)==1) or (self.string[1] not in ['-','+'])):
      #    # one of the four bases
      #    self.types[self.string[0]] += 1
      #    self.string = self.string[1:]
          
      else:
          # unrecognized character
          # or a read that reports a substitition followed by an insertion/deletion
          self.types['X'].append(self.string[0])
          self.string = self.string[1:]
          self.types['Reads'] += 1
          self.types['NotSNP'] += 1
          
    return
  
  def __repr__(self):
    types = self.types
    return '\t'.join( list(map(str,[types['A'], types['C'], types['G'],types['T'], 
                                    types['*'], types['NotSNP'],types['SNP_Reads'],types['Reads']])) + 
                      list(map(','.join, [types['-'],types['+'],types['X']])))
    

def main():
  
  sys.stderr.write("NOTICE: the parser does not count the match nt allele in cases of adjacent insertion or deletion.\n" + \
                     " This parser is valid only when insertions / deletions are <=9bp long.\n")
  
  sys.stdout.write("Chr\tPosition\tref\tcov\tA\tC\tG\tT\tstar\tNotSNP\tSNP_Reads\tReads\tinsertions\tdeletions\tunknown\n")
  for line in sys.stdin:  
    #sys.stderr.write("DEBUG line: %s\n" % line.strip())
    toks = line.strip('\n').split('\t')
    ref = toks[2].upper()
    cov = toks[3]
    sys.stdout.write('\t'.join([toks[0], toks[1],ref, cov]) + '\t' + \
        parseString(ref, toks[4]).__repr__() + '\n')

if __name__ == '__main__':
  main()