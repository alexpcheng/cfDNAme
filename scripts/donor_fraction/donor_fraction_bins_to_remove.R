library(data.table)

args = commandArgs(trailingOnly = TRUE)

file_path = args[1]
out_file = args[2]

max_kept_cov = args[3]

file.list<-list.files(path = file_path, pattern = "depth_binned_cov.txt", recursive = TRUE)

fread_path <- function(file, path){
  return(fread(paste0(path, file), header=FALSE))
}

df<-rbindlist(lapply(X = file.list, FUN=fread_path, file_path))

aggregate_coverage<-aggregate(.~V1+V2+V3, df, sum)
colnames(aggregate_coverage)<-c("chrom", "binstart", "binend", "summed_coverage")
aggregate_coverage <- aggregate_coverage[aggregate_coverage$summed_coverage >= max_kept_cov, ]
aggregate_coverage <- aggregate_coverage[order(aggregate_coverage$chrom, aggregate_coverage$binstart), ]
fwrite(x = aggregate_coverage[, c("chrom", "binstart"), drop=FALSE], 
       file=out_file, 
       quote = FALSE, 
       row.names = FALSE, 
       col.names = TRUE,
       sep = '\t')
