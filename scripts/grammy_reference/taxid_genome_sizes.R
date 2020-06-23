library(data.table)

args=commandArgs(trailingOnly = TRUE)
fai_file = args[[1]]
gi_taxid_file = args[[2]]

fai = fread(fai_file, col.names = c('gi', 'genome_len', 'a', 'b', 'c'), sep = '\t')
gi_taxid = fread(gi_taxid_file, col.names = c('gi', 'taxid'))
output = args[[3]]


fai$gi <- sapply(X = fai$gi, FUN=function(x) strsplit(x = x, split = '[|]')[[1]][2])
fai$gi <- as.numeric(fai$gi)

gi_taxid$gi <- as.numeric(gi_taxid$gi)

taxid_len <- merge(fai[, c('gi', 'genome_len')],
                   gi_taxid, by = 'gi')

taxid_len <- taxid_len[, c('taxid', 'genome_len')]
taxid_len <- aggregate(.~taxid, taxid_len, max) #used to be mean
taxid_len$genome_len <- round(taxid_len$genome_len)
fwrite(x = taxid_len, file=output, quote = FALSE, col.names = TRUE, sep='\t')
