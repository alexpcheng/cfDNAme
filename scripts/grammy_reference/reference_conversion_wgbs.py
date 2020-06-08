#!/usr/bin/env python
import sys

def parse_arguments():
    input_fasta = sys.argv[1]
    CT = sys.argv[2]
    GA = sys.argv[3]
    read_ids = sys.argv[4]
    return(input_fasta, CT, GA, read_ids)

def convert_database(input_fasta, CT, GA, read_ids):
    with open(input_fasta) as f, open(CT, 'w') as ct, open(GA, 'w') as ga, open(read_ids, 'w') as ri:
        gis = []
        for line in f:
            if line[0] == '>':
                ct.write(line)
                ga.write(line)
                gi=line.split('|')[1]
                gis.append(gi)
            else:
                ct.write(line.upper().replace('C', 'T'))
                ga.write(line.upper().replace('G', 'A'))
        ri.write('\n'.join(gis))

def main():
    input_fasta, CT, GA, read_ids = parse_arguments()
    convert_database(input_fasta, CT, GA, read_ids)

if __name__ == '__main__':
    main()
