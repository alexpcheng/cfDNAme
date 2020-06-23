library(data.table)

get_eff_len2 <- function(a, b, c, d, e){
  a = as.numeric(a)
  b = as.numeric(b)
  c = as.numeric(c)
  d = as.numeric(d)
  genome_len = as.numeric(e)
  lengths = sort(c(a,b,c,d))
  len1 = lengths[4]-lengths[1]
  len2 = genome_len - lengths[3] + lengths[2]
  eff_len = min(len1, len2)
  return(eff_len)
}

vv <- Vectorize(get_eff_len2, vectorize.args=c('a', 'b', 'c', 'd', 'e'))


args = commandArgs(trailingOnly=TRUE)
R1_file = args[[1]]
R2_file = args[[2]]
outfile = args[[3]]
taxid_file = args[[4]]

taxid_lengths = fread(taxid_file)
taxid_lengths$taxid<- as.character(taxid_lengths$taxid)

R1 = fread(R1_file, col.names = c('qseqid', 'sseqid_R1', 'pident_R1', 'length_R1', 'mismatch_R1', 'gapopen_R1', 'qstart_R1', 'qend_R1', 'sstart_R1', 'send_R1', 'evalue_R1', 'bitscore_R1', 'qlen_R1', 'strand', 'taxid'))
R2 = fread(R2_file, col.names = c('qseqid', 'sseqid_R2', 'pident_R2', 'length_R2', 'mismatch_R2', 'gapopen_R2', 'qstart_R2', 'qend_R2', 'sstart_R2', 'send_R2', 'evalue_R2', 'bitscore_R2', 'qlen_R2', 'strand', 'taxid'))

R1 = R1[R1$length_R1/R1$qlen_R1 >=0.9, ]
R2 = R2[R2$length_R2/R2$qlen_R2 >=0.9, ]

R1 = R1[R1$taxid != "GI_NOT_FOUND", ]
R2 = R2[R2$taxid != "GI_NOT_FOUND", ]

R1$taxid = as.character(R1$taxid)
R2$taxid = as.character(R2$taxid)

R1$qseqid = gsub("_.*", "", R1$qseqid)
R2$qseqid = gsub("_.*", "", R2$qseqid)
R1 = R1[(R1$qseqid %in% R2$qseqid), ]
R2 = R2[(R2$qseqid %in% R1$qseqid), ]

#R1 = R1[(R1$evalue_R1<10**-20 & R1$bitscore_R1>100), ]
#R2 = R2[(R2$evalue_R2<10**-20 & R2$bitscore_R2>100), ]

read_ids <- unique(R1$qseqid)

over_occuring_rows <- merge(data.frame(table(R1$qseqid)), 
                            data.frame(table(R2$qseqid)), 
                            by='Var1')
over_occuring_reads <- as.character(over_occuring_rows$Var1[over_occuring_rows$Freq.x >10000 | over_occuring_rows$Freq.y >10000])

R1 <- R1[!(R1$qseqid %in% over_occuring_reads), ]
R2 <- R2[!(R2$qseqid %in% over_occuring_reads), ]

paired = merge(R1, R2, by=c('qseqid', 'strand', 'taxid'), allow.cartesian = TRUE)

paired = merge(paired, taxid_lengths, by='taxid')

paired$effective_length <- with(paired, vv(sstart_R1, send_R1, sstart_R2, send_R2, genome_len))

paired = paired[paired$effective_length<1000, ]

paired = paired[order(paired$qseqid, paired$effective_length), ] #duplicate function will keep the first occurence, which by this definition is
# also the shortest
to_keep = !duplicated(paired$qseqid)

paired = paired[to_keep, ]

fwrite(x = paired, file = outfile, quote = FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
