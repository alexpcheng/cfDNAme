library(data.table)
args = commandArgs(trailingOnly = TRUE)
input = args[[1]]
output = args[[2]]
sam = args[[3]]
df <- fread(input, header=FALSE, data.table = FALSE)
df <- data.frame(t(df))
df <- df[df$X2>0, ]
df$SAM <- sam
df <- df[, c('SAM', 'X1', 'X2', 'X3')]
fwrite(x = df, file=output, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')