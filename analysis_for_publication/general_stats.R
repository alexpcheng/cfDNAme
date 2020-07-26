#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize workspace -----------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)

workdir <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/'
setwd(workdir)

master_table_path <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/master_table_megakaryocytes.csv'
figure_path <- '/workdir/apc88/WGBS_pipeline/Figures/'

# Read in master table -----------------------------------------------------------------------
df <- fread(master_table_path, na.strings = "NA")

# STATISTICS ---------------------------------------------------------------------------------
# COVID-19+ plasma samples
table(df$sample_class)
table(df$patient_id)

# COVID-19+ patients
length(unique(df[df$sample_class=="SARS-CoV2+", ]$patient_id))

# UCSF covid-19+ plasma samples
nrow(df[(df$cohort == "high_frequency_collection") & (df$sample_class == "SARS-CoV2+"), ])

# MCGILL covid 19+ plasma samples
nrow(df[(df$cohort == "randomized_control") & (df$sample_class == "SARS-CoV2+"), ])

# MCGILL S-O-C
length(unique(df[(df$random_trial_treatment == "standard_of_care"), ]$patient_id))
length(unique(df[(df$random_trial_treatment == "lopinavir/ritonavir"), ]$patient_id))

# MCGILL patient hospitalization
table(unique(df[df$cohort=="randomized_control", c("patient_id", "hospitalization_status_as_of_june_26_2020")])$hospitalization_status_as_of_june_26_2020)

table((df[df$cohort=="randomized_control", c("patient_id", "hospitalization_status_as_of_june_26_2020")])$hospitalization_status_as_of_june_26_2020)


# other flu patients
length(unique(df[df$sample_class=="other_RNA_virus", ]$patient_id))

# total PE reads
mean(df$total_number_of_pe_reads)/10**6
sd(df$total_number_of_pe_reads)/10**6

# coverage
mean(df$hg19_cov)
sd(df$hg19_cov)

# bisulfite
mean(df$bisulfite_conversion)
sd(df$bisulfite_conversion)

