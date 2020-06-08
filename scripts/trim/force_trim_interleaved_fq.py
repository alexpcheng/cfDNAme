import sys

R1_out = sys.argv[1]
R2_out = sys.argv[2]
min_len = 10 #matches bbduk
with open(R1_out, 'w') as r1, open(R2_out, 'w') as r2:
    dat = []
    for line in sys.stdin:
        dat.append(line)
        if len(dat) == 8: #reads the the two interleaved reads
            rid_1 = dat[0]
            seq_1 = dat[1]
            plu_1 = dat[2]
            qua_1 = dat[3]

            rid_2 = dat[4]
            seq_2 = dat[5][15:]
            plu_2 = dat[6]
            qua_2 = dat[7][15:]

            if ( len(seq_2.rstrip('\n')) > min_len and len(seq_1.rstrip('\n')) > min_len):
                r1.write(rid_1 + seq_1 + plu_1 + qua_1)
                r2.write(rid_2 + seq_2 + plu_2 + qua_2)
            dat = []
