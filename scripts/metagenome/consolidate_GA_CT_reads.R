library(data.table)

args= commandArgs(trailingOnly = TRUE)

GA_file = args[[1]] 
CT_file = args[[2]]
tblatPE_file = args[[3]]
tblat1_file = args[[4]]

GA = fread(GA_file, col.names = c("taxid", "qseqid", "strand", "sseqid_R1", "pident_R1","length_R1", "mismatch_R1", "gapopen_R1", "qstart_R1",
                                  "qend_R1", "sstart_R1", "send_R1", "evalue_R1", "bitscore_R1", "qlen_R1", "sseqid_R2", "pident_R2", "length_R2", 
                                  "mismatch_R2", "gapopen_R2", "qstart_R2", "qend_R2", "sstart_R2", "send_R2", "evalue_R2", "bitscore_R2", "qlen_R2",
                                  "genome_len", "effective_length"))
CT = fread(CT_file, col.names = c("taxid", "qseqid", "strand", "sseqid_R1", "pident_R1","length_R1", "mismatch_R1", "gapopen_R1", "qstart_R1",
                                  "qend_R1", "sstart_R1", "send_R1", "evalue_R1", "bitscore_R1", "qlen_R1", "sseqid_R2", "pident_R2", "length_R2", 
                                  "mismatch_R2", "gapopen_R2", "qstart_R2", "qend_R2", "sstart_R2", "send_R2", "evalue_R2", "bitscore_R2", "qlen_R2",
                                  "genome_len", "effective_length"))

tblatPE = rbind(GA, CT)
tblatPE = tblatPE[order(tblatPE$qseqid, tblatPE$genome_len), ]

to_keep = !duplicated(tblatPE$qseqid)

tblatPE = tblatPE[to_keep, ]
tblat1 = tblatPE[, c('qseqid', 'sseqid_R1', 'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', 'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', 'strand', 'taxid')]
tblat1$qseqid <- paste0(tblat1$qseqid, '-1')
fwrite(x = tblatPE, file = tblatPE_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
fwrite(x = tblat1, file = tblat1_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)