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
# Functions for tissue of origin measurement
group_by_celltype <- function(df){
  tmpA <- data.frame(t(df))
  tmpA$tissue <- gsub('_.*', '', colnames(df))
  new_df <- aggregate(.~tissue, tmpA, mean)
  new_names <- new_df$tissue
  new_df$tissue<-NULL
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- new_names
  new_df <- as.matrix(new_df)
  return(new_df)
}

# mixing parameters function
mp_function <- function(ref_and_sam, other, sum_to_one, group_tissues, sample_name){
  A <- ref_and_sam[, !(c('chr', 'start', 'end', 'pmeth'))]
  b <- ref_and_sam$pmeth
  
  if (other){
    A$other <- 1
  }
  
  if (group_tissues){
    A <- group_by_celltype(A)
  }else{
    A <- as.matrix(A)
  }
  
  mod <- nnls(A, b)
  percentages <- data.frame(mod$X)
  
  colnames(percentages) <-c('prediction')
  percentages$prediction <- as.numeric(as.character(percentages$prediction))
  
  if (sum_to_one){
    percentages$prediction <- percentages$prediction/sum(percentages$prediction)
  }
  percentages <- data.frame(t(percentages))
  colnames(percentages) <- colnames(A)
  percentages$sample <- sample_name
  return(percentages)
  
}

# Command line arguments -------------------------------------------------------------------------------------

# ARGUMENTS MUST BE SPECIFIED IN THE FOLLOWING ORDER
# 1- sample file path
# 2- reference file path
# 3- output file path
# 4- lookup table (reference methylomes)
# 5- sample name
# 6- sum to one: do we want to sum to one? TRUE or FALSE
# 7- other: do we want the possibility of an error absorbing term? TRUE or FALSE
# 8- group by cell type: do we average the references by group first? TRUE or FALSE
# removals: what tissues in the reference to you want to remove before deconvoluting? no spaces, front slash separated
#               ex: bcell_1/t_cell3

args <- commandArgs(trailingOnly = TRUE)

sam <- args[[1]]
refs <- args[[2]]
out <- args[[3]]
lookup <- args[[4]]
sample_name <- args[[5]]
sum_to_one <- args[[6]]
other <- args[[7]]
group_by_celltype <- args[[8]]


if (!(sum_to_one == "TRUE" | sum_to_one == "FALSE")){
  stop('Your sum_to_one parameter is wrong. TRUE or FALSE with caps.')
}else{
  sum_to_one <- as.logical(sum_to_one)
}
if (!(other == "TRUE" | other == "FALSE")){
  stop('Your other parameter is wrong. TRUE or FALSE with caps.')
}else{
  other <- as.logical(other)
}
if (!(group_by_celltype == "TRUE" | group_by_celltype == "FALSE")){
  stop('Your group_by_celltype parameter is wrong. TRUE or FALSE with caps.')
}else{
  group_by_celltype <- as.logical(group_by_celltype)
}

if ((length(args)) ==9){
  removals <- args[[9]]
  removals <- gsub('/', '|', removals)
}else{
  removals <- 'thiswillneverbeatissuename' 
}

# Read files and format column names -----------------------------------------------------------------------
sam.df <- fread(sam)
colnames(sam.df)<- c('chr', 'start', 'end', 'pmeth')
sam.df$pmeth <- sam.df$pmeth/100

refs.df <- fread(refs)
coln<-fread(lookup, header=T)
colnames(refs.df)<-c("chr", "start", "end", coln$tissue_name)

refs.df <- data.table(data.frame(refs.df)[, !grepl(removals, colnames(refs.df))])

# Initialize variables ----------------------------------------------------------------------------------------
num.tissues<-ncol(refs.df)-3 # 3 columns are for chromosome, start and end

# Table preparation ------------------------------------------------------------------------------------------
ref_and_sam <- merge(refs.df, sam.df, by=c('chr', 'start', 'end'))
ref_and_sam <- ref_and_sam[complete.cases(ref_and_sam)]

# Measuring tissue of orign ----------------------------------------------------------------------------------
tissues_of_origin <- mp_function(ref_and_sam = ref_and_sam,
                                 other=other,
                                 sum_to_one = sum_to_one,
                                 group_tissues = group_by_celltype,
                                 sample_name=sample_name)


#tissues_of_origin
fwrite(tissues_of_origin, file = out, quote=FALSE, sep="\t")

