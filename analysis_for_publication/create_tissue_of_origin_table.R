#!/usr/bin/env RScript
# Title: UCSF_TOO.R
# Authors: Alexandre Pellan Cheng (apc88@cornell.edu)
# Brief description: Analysis of tissues of origin from UCSF cohort (high frequency collection)

# Initialize workspace --------------------------------------------------------------------------
rm(list = ls())

library(data.table)

# Get QC metrics for filtering -----------------------------------------------------------------
cov <- rbindlist(lapply(X=list.files(path = '/workdir/apc88/WGBS_pipeline/sample_output/stats/',
                                     pattern = 'mSARS|14',
                                     full.names = TRUE),
                        FUN=function(x){df<-fread(x); df$sample<-strsplit(x = x, split='/')[[1]][8]; return(df)}))
cov$sample <- gsub('_mapping_stats.txt', '', cov$sample)
cov <- cov[order(cov$sample, cov$NUM_READS), ]

tot_ab <- rbindlist(lapply(X=list.files(path = '/workdir/apc88/WGBS_pipeline/sample_output/total_abundance/',
                                        pattern = 'mSARS|14',
                                        full.names = TRUE),
                           FUN = fread))

df <- merge(cov, tot_ab)
df$sample_id <- gsub('_S.*', '', df$sample)
df$sample_id <- gsub("m$", '', df$sample_id)

pure_pass <- df[df$DEPTH > 0.2 & df$spikeins >=10, c('sample', 'sample_id', 'DEPTH', 'spikeins', 'human_aligned')]

tissues <- rbindlist(lapply(X= list.files(path = '/workdir/apc88/WGBS_pipeline/sample_output/MEGAKA_Refs/classic_tissues_of_origin_trim/',
                                          pattern = '.tsv', full.names = TRUE),
                            FUN = fread))

tissues <- tissues[grepl("MCGILL|SWIFTHLY", tissues$sample) | tissues$sample %in% pure_pass$sample,]

tissues$other <- NULL

tissues.melt <- melt(tissues, id.vars = 'sample')
tissues.melt$variable <- gsub("_.*", "", tissues.melt$variable)
tissues.melt <- aggregate(.~variable+sample, tissues.melt, sum)

tissues.melt$variable <- factor(tissues.melt$variable, levels = c('macrophage', 'monocyte', 'dendritic',
                                                                  'neutrophil', 'eosinophil', 
                                                                  'tcell', 'nkcell', 'bcell',
                                                                  'erythroblast', 'hsc', 'progenitor',
                                                                  'spleen', 'kidney', 'heart', 'skin', 'lung',
                                                                  'liver', 
                                                                  'pancreas',
                                                                  'colon', 'megakaryocyte'))
df <- dcast(data = tissues.melt, formula = sample~variable)

df <- df[!(df$sample %in% c("MCGILL1", "MCGILL23", "MCGILL26", "MCGILL35", "MCGILL36")), ]
df <- df[!grepl("DUP", df$sample), ]
df$sample <- gsub("CAT", "", df$sample)
df <- df[order(df$sample), ]
