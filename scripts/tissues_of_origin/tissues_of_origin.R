#!/usr/bin/env RScript
# Title: Tissue of Origin Measurements
# Authors: Alexandre Pellan Cheng
# Brief description: Creates a text file with the tissue of origin measurements (aka mixing parameters (mp)) for
# sample. References are loaded and have already assumed to be prepped (via references.500.rds)

# setwdInitialize -------------------------------------------------------------------------------------------------
rm(list=ls())
# Load libraries ---------------------------------------------------------------------------------------------
library(Matrix)
library(matrixcalc)
library(stats)
library(data.table)
library(limSolve)
library(parallel)

# Source custom function ------------------------------------------------------------------------------------- 
source("scripts/tissues_of_origin/tissues_of_origin_functions.R")

# Command line arguments -------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
print('waaaa')
sam <- args[[1]]
refs <- args[[2]]
out <- args[[3]]
lookup <- args[[4]]
sample_name <- args[[5]]
# Read files and format column names -----------------------------------------------------------------------
sam.df <- fread(sam)
colnames(sam.df)<- c('chr', 'start', 'end', 'pmeth')
sam.df$pmeth <- sam.df$pmeth/100

refs.df <- fread(refs)
coln<-fread(lookup, header=T)
colnames(refs.df)<-c("chr", "start", "end", coln$tissue_name)

##
#n = 30
#pcref <- prcomp(t(refs.df[, 4:ncol(refs.df)]))
#Xhat = pcref$x[, 1:n] %*% t(pcref$rotation[, 1:n])
#Xhat = scale(Xhat, center=-colMeans(t(refs.df[, 4:ncol(refs.df)])), scale=FALSE)
#r <- data.frame(refs.df[, c(1:3)], t(Xhat))
#refs.df <- r
#colnames(refs.df)<-c("chr", "start", "end", coln$V3)
##

# Initialize variables ----------------------------------------------------------------------------------------
num.tissues<-ncol(refs.df)-3 # 3 columns are for chromosome, start and end

# Table preparation ------------------------------------------------------------------------------------------
ref_and_sam <- merge(refs.df, sam.df, by=c('chr', 'start', 'end'))
ref_and_sam <- ref_and_sam[complete.cases(ref_and_sam)]

# Measuring tissue of orign ----------------------------------------------------------------------------------
tissues_of_origin <- mp_function(ref_and_sam, other=TRUE, nt=num.tissues, sum_to_one = TRUE, group_tissues = TRUE, sample_name=sample_name)
#tissues_of_origin
fwrite(tissues_of_origin, file = out, quote=FALSE, sep="\t")

