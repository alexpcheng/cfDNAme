#!/usr/bin/env python3
import sys
def parse_args():
    wig = sys.argv[1]
    out = sys.argv[2]
    chr = sys.argv[3:]
    return(wig, out, chr)

def filter(wig, out, chr):
    with open(wig) as f, open(out, 'w') as w:
        for line in f:
            if 'fixed' in line:
                chromosome = line.split(' ')[1].split('=')[1]
                if chromosome in chr:
                    to_write=True
                else:
                    to_write=False
            if to_write:
                w.write(line)

def main():
    wig, out, chr = parse_args()
    filter(wig, out, chr)

if __name__ == '__main__':
    main()
