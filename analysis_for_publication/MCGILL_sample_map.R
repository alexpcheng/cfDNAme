#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize workspace -----------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)
library(stringr)
source('~/theme_alex.R')

workdir <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/'
setwd(workdir)

master_table_path <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/master_table_megakaryocytes.csv'
figure_path <- '/workdir/apc88/WGBS_pipeline/Figures/'

# Read in master table -----------------------------------------------------------------------
df <- fread(master_table_path, na.strings = "NA")

# Longitudinal analysis of COVID19 -----------------------------------------------------------

MCGILL_df <- df[df$cohort == "randomized_control", ]
MCGILL_df$rename <- as.character(as.numeric(gsub("MED", "", MCGILL_df$patient_id))/10) #formats for axis sizes to match numerical data frames
MCGILL_df$hospitalization_status_as_of_june_26_2020 <-
  factor(MCGILL_df$hospitalization_status_as_of_june_26_2020, levels = c("discharged", "hospitalized", "deceased"))

agg1 <- MCGILL_df
agg1$c <- 1

agg1 <- aggregate(.~patient_id, agg1[, c('patient_id', 'c')], sum)

agg2 <- aggregate(.~patient_id, MCGILL_df[, c('patient_id', 'days_past_trial_enrolment')], min)

MCGILL_df <- merge(MCGILL_df, agg1, by = 'patient_id')
MCGILL_df <- merge(MCGILL_df, agg2, by = 'patient_id', all.x = TRUE)

MCGILL_df <- MCGILL_df[order(MCGILL_df$hospitalization_status_as_of_june_26_2020,
                             MCGILL_df$c,
                             MCGILL_df$days_past_trial_enrolment.y, decreasing = TRUE), ]

MCGILL_df$rename <- factor(MCGILL_df$rename, levels = MCGILL_df$rename[!duplicated(MCGILL_df$rename)])

MCGILL_df$ord <- NA
for (i in 1:nrow(MCGILL_df)){
  if (MCGILL_df$random_trial_timepoint[i] == "D0"){
    MCGILL_df$ord[i]<- MCGILL_df$ordinal_D0[i]
  }
  if (MCGILL_df$random_trial_timepoint[i] == "D1"){
    MCGILL_df$ord[i]<- MCGILL_df$ordinal_D0[i]
  }
  if (MCGILL_df$random_trial_timepoint[i] == "D5"){
    MCGILL_df$ord[i]<- MCGILL_df$ordinal_D5[i]
  }
  if (MCGILL_df$random_trial_timepoint[i] == "D15"){
    MCGILL_df$ord[i]<- MCGILL_df$ordinal_D15[i]
  }
}

# FIGURES ------------------------------------------------------------------------------------
if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_sample_map.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data=MCGILL_df, aes(x=days_past_trial_enrolment.x, y= rename, color = hospitalization_status_as_of_june_26_2020))+
  geom_line(aes(group = patient_id), linetype = 'dashed')+
  geom_point(aes(shape = random_trial_treatment))+
  scale_color_manual(values = c('purple', 'deepskyblue1', 'red'))+
  scale_shape_manual(values=c(1,2))+
  xlab('Days post trial entry')+
  ylab('Patient')+theme_alex()+
  #theme(axis.text.y=element_blank())+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_sample_map2.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data=MCGILL_df, aes(x=days_past_trial_enrolment.x, y= rename))+
  geom_line(aes(group = patient_id), linetype = 'dashed')+
  geom_point(aes(shape = hospitalization_status_as_of_june_26_2020, 
                 fill = as.numeric(as.character(ord))))+
  #scale_color_manual(values = c('purple', 'deepskyblue1', 'red'))+
  scale_shape_manual(values=c(21,22, 24))+
  scale_fill_viridis_c()+
  xlab('Days post trial entry')+
  ylab('Patient')+
  xlim(c(-1,21))+
  theme_alex()+
  #theme(axis.text.y=element_blank())+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_sample_map2_LEG.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data=MCGILL_df, aes(x=days_past_trial_enrolment.x, y= rename))+
  geom_line(aes(group = patient_id), linetype = 'dashed')+
  geom_point(aes(shape = hospitalization_status_as_of_june_26_2020, 
                           fill = as.numeric(as.character(ord))))+
  #scale_color_manual(values = c('purple', 'deepskyblue1', 'red'))+
  scale_shape_manual(values=c(21,22, 24))+
  scale_fill_viridis_c()+
  xlab('Days post trial entry')+
  ylab('Patient')+
  xlim(c(-1,21))+
  #theme(axis.text.y=element_blank())+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}