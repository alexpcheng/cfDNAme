args=commandArgs(trailingOnly = TRUE)

grammy.file <- args[[1]]
blast.file <- args[[2]]
stats.file <- args[[3]]
REF <- args[[4]]
output.file <- args[[5]]
grammy.tab <- read.table(grammy.file, header = FALSE, fill = TRUE)
colnames(grammy.tab) <- c("SAMPLE", "Taxid", "GrAb", "GrEr")
#
blast <- read.table(blast.file, header = FALSE, fill = TRUE)
total.blast <- nrow(blast)
#
align.stats <- read.table(stats.file, header = TRUE, fill = TRUE)
hg.coverage <- align.stats$DEPTH[1]
#
grammy.LUT <- read.table(REF, header = TRUE, fill = TRUE)

grammy.tab.info <- merge(grammy.tab, grammy.LUT, by = "Taxid")
# weighted genome size
grammy.tab.info$hgcoverage <- hg.coverage
grammy.tab.info$WeightedGenome <- sum(grammy.tab.info$Length * grammy.tab.info$GrAb)
grammy.tab.info$AdjustedBlast <- total.blast*(grammy.tab.info$Length*grammy.tab.info$GrAb/grammy.tab.info$WeightedGenome)
grammy.tab.info$Coverage <- 75*grammy.tab.info$AdjustedBlast/grammy.tab.info$Length
grammy.tab.info$RelCoverage <- 2*75*grammy.tab.info$AdjustedBlast/grammy.tab.info$Length/hg.coverage
#

write.table(grammy.tab.info, output.file,sep ="\t", row.names = FALSE)
