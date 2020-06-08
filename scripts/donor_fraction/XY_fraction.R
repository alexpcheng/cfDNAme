### Burnham, De Vlaminck: 2017
# donorfraction.autosomalandY.R

#Slight modifications from original script done by Alexandre Pellan Cheng
#Adapted for BS treated samples
### libraries
suppressMessages(library(HMMcopy))
source('scripts/donor_fraction/HMM_modif.R')
library(data.table)
### Variables
args <- commandArgs(trailingOnly = T) # command line arguments

sam = as.character(args[1])

rfile = args[2]
mfile = args[3]
gfile = args[4]
outfile = args[5]
removals <- args[6]
####
removals <- fread(removals)
removals$remove<-"Yes"
####

# Function -----------------------------------------------------------------------------------------------------------
get_sex_frac <- function(sam, rfile, mfile, gfile, thr=0.8, stdv_multi = 6){
  normal_reads <- wigsToRangedData(readfile = rfile, gcfile = gfile,mapfile = mfile, verbose=F)
  normal_copy <- correctReadcount2(normal_reads,verbose = F, samplesize = Inf)
  #change to easy to read data.frame format
  
  #### TEST OLD VERSION #####
  df <- cbind(data.frame(normal_copy@ranges),normal_copy$cor.map, normal_copy$map) #, data.frame(normal_copy@values))
  #df <- df[df$ideal, ]
  
  df <- df[,c('group_name', 'start', 'normal_copy$cor.map', 'normal_copy$map')]
  colnames(df)<-c("chrom","binstart","value","mapa")
  #df[is.na(df)] <- 0
  df<-df[!is.na(df$value), ]
  #determine cutoffs to eliminate overrepresented bins
  ythresh <- c(mean(df[df$chrom=="chrY"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom=="chrY"&!is.na(df$value),]$value),
               mean(df[df$chrom=="chrY"&!is.na(df$value),]$value) - stdv_multi*sd(df[df$chrom=="chrY"&!is.na(df$value),]$value))
  athresh <- c(mean(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value),
               mean(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value) - stdv_multi*sd(df[df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value))
  xthresh <- c(mean(df[df$chrom=="chrX"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom=="chrX"&!is.na(df$value),]$value),
               mean(df[df$chrom=="chrX"&!is.na(df$value),]$value) - stdv_multi*sd(df[df$chrom=="chrX"&!is.na(df$value),]$value))
  
  chr1thresh<-c(mean(df[df$chrom=="chr1"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom=="chr1"&!is.na(df$value),]$value),
                mean(df[df$chrom=="chr1"&!is.na(df$value),]$value) - stdv_multi*sd(df[df$chrom=="chr1"&!is.na(df$value),]$value))
  
  a_no_chr1thresh <- c(mean(df[df$chrom!="chr1" & df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value) + stdv_multi*sd(df[df$chrom!="chr1" & df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value),
                       mean(df[df$chrom!="chr1" & df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value) - stdv_multi*sd(df[df$chrom!="chr1" & df$chrom!="chrX"&df$chrom!="chrY"&!is.na(df$value),]$value))
  
  #### ADDED. SHOULD BE THE SAME AS COMMENTED LINE
  ysa <- df[df$chrom == "chrY" & df$mapa > thr & df$value <ythresh[1] &df$value>ythresh[2], ]
  ysa <- merge(as.data.table(ysa), removals, by=c("chrom", "binstart"), all.x=TRUE)
  ysa <- data.frame(ysa[is.na(ysa$remove), c("chrom", "binstart", "value", "mapa")])
  ###
  
  asa <- df[df$chrom!="chrX" & df$chrom!="chrY" & df$mapa > thr & df$value <athresh[1] & df$value > athresh[2],]
  #### ADDED NEW FEATURE. Y-DERIVED CFDNA THAT MAPS TO AUTOSOMALS
  asa <- merge(as.data.table(asa), removals, by = c("chrom", "binstart"), all.x=TRUE)
  asa <- data.frame(asa[is.na(asa$remove), c("chrom", "binstart", "value", "mapa")])
  ####

  xsa <- df[df$chrom == "chrX" & df$mapa > thr & df$value <xthresh[1] & df$value >xthresh[2], ]
  xsa <- merge(as.data.table(xsa), removals, by=c("chrom", "binstart"), all.x=TRUE)
  xsa <- data.frame(xsa[is.na(xsa$remove), c("chrom", "binstart", "value", "mapa")])
  
  chr1sa <- df[df$chrom=="chr1" & df$mapa > thr & df$value < chr1thresh[1] & df$value > chr1thresh[2],]
  chr1sa <- merge(as.data.table(chr1sa), removals, by=c("chrom", "binstart"), all.x=TRUE)
  chr1sa <- data.frame(chr1sa[is.na(chr1sa$remove), c("chrom", "binstart", "value", "mapa")])
  
  no_chr1sa <- df[df$chrom!="chr1" & df$chrom!="chrX" & df$chrom!="chrY" & df$mapa > thr & df$value <a_no_chr1thresh[1] & df$value > a_no_chr1thresh[2],]
  no_chr1sa <- merge(as.data.table(no_chr1sa), removals, by=c("chrom", "binstart"), all.x=TRUE)
  no_chr1sa <- data.frame(no_chr1sa[is.na(no_chr1sa$remove), c("chrom", "binstart", "value", "mapa")])
  #calculate donor fraction based on the relative coverages of chrY and the average autosomal coverage.
  
  X <- 100*mean(xsa$value)/mean(asa$value)
  Y <- 100*mean(ysa$value)/mean(asa$value)
  chr1 <- 100*mean(chr1sa$value)/mean(no_chr1sa$value)

  df <- data.frame(sam, X, Y, chr1)
  return(df)
}

## Code execution #################################################

result = get_sex_frac(sam = sam, rfile=rfile, gfile=gfile, mfile=mfile)

fwrite(x=result, file=outfile, quote=FALSE, sep='\t', na = '0')
