#!/usr/bin/env RScript
# Title: Histogram Fragment Lengths
# Authors: Alexandre Pellan Cheng, Iwijn De Vlaminck
# Brief description: Creates fragment length profiles for paired end sequencing data

# Initialize -------------------------------------------------------------------------------------------------
rm(list=ls())
source('~/theme_alex.R')
# Load libraries ------------------------------------------------------------------------
require(ggplot2)
library(data.table)
# Initialize variables ------------------------------------------------------------------
save_eps <- T

args = commandArgs(trailingOnly=TRUE)
LengthsFile <- args[1]
Fig <- args[2]

# Read and format length files ----------------------------------------------------------
lengths <-  data.frame(fread(LengthsFile, header = FALSE, col.names = c('count', 'fragment_length')))
lengths$frac <- lengths$count/sum(lengths$count)

# Plot and save data --------------------------------------------------------------------
if (save_eps) { pdf(file=Fig, width=6.5/2.5,height=5/2.5, paper="special",bg="white",
                    fonts="Helvetica", colormodel="cmyk", pointsize = 8)}

ggplot(lengths, aes(y=frac, x=fragment_length)) +
  geom_line(colour = "red")+ 
  xlim(0,max(lengths$fragment_length)) + 
  xlab('Fragment length (bp)')+ylab('Fraction of reads')+
  theme_alex()

if (save_eps) {dev.off()}
