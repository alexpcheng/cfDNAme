#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

## Loading packages
library(ggplot2)

## make nice plots
theme_cfg <-
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),

    axis.line = element_line(colour = "black"),

    legend.key = element_blank()
  )

## Read input file
data <- read.table(file=args[1], header=F, col.names=c('chr','start','end','q_val','diff','CpG','group1','group2'))
data$length <- data$end-data$start

## Plot statistics
pdf(args[2])
#difference histogram
ggplot(data, aes(x=diff)) + geom_histogram(binwidth=0.02, fill='darkgrey', color='black') + xlab("mean methylation difference") + theme_cfg + scale_x_continuous(limits=c(-1,1))
#length distribution CpGs
ggplot(data, aes(x=length)) + geom_line(stat="density", size=1) + xlab("DMR length [nt]") + theme_cfg
#length distribution nt
ggplot(data, aes(x=CpG)) + geom_line(stat="density", size=1) + xlab("DMR length [CpG]") + theme_cfg
#q_val vs difference
ggplot(data, aes(x=diff, y=q_val)) + geom_point(alpha=.5) + scale_y_log10() + xlab("mean methylation difference") + ylab("q-value") + theme_cfg
#mean1 vs mean2
ggplot(data, aes(x=group1, y=group2)) + geom_point(alpha=.5) + coord_fixed() + xlab("mean group 1") + ylab("mean group 2") + theme_cfg
#nt vs CpGs
ggplot(data, aes(x=length, y=CpG)) + geom_point(alpha=.5) + xlab("DMR length [nt]") + ylab("DMR length [CpGs]") + theme_cfg
dev.off()

