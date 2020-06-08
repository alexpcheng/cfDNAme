library(Rsamtools)
library(data.table)

args = commandArgs(trailingOnly = TRUE)


full_stats = fread(args[[2]])

num_pe_spikeins = countBam(args[[1]])$records/2

df <- data.frame(args[[4]], num_pe_spikeins, full_stats$ALIGNED)

colnames(df)<-c('sample', 'spikeins', 'human_aligned')

fwrite(df, args[[3]], quote = FALSE, sep='\t')
