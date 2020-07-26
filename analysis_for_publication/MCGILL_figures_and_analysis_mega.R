#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize workspace -----------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)
source('~/theme_alex.R')
workdir <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/'
setwd(workdir)

master_table_path <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/master_table_megakaryocytes.csv'
figure_path <- '/workdir/apc88/WGBS_pipeline/Figures/'

# Read in master table -----------------------------------------------------------------------
df <- fread(master_table_path, na.strings = "NA", data.table = FALSE)
df$Notes<-NULL

tissues <- c('macrophage', 'monocyte', 'dendritic','neutrophil', 'eosinophil', 'Tcell', 'Nkcell', 'Bcell',
             'erythroblast', 'hsc', 'progenitor','spleen', 'kidney', 'heart', 'skin', 'lung','liver','pancreas',
             'colon', 'megakaryocyte')
# Longitudinal analysis of COVID19 -----------------------------------------------------------

MCGILL_df <- df[df$cohort == "randomized_control" | df$cohort=="healthy_control", ]
MCGILL_df[is.na(MCGILL_df)]<-"HLY"
MCGILL_df$random_trial_timepoint[MCGILL_df$random_trial_timepoint=="D1"]<-"D0"
MCGILL_df$random_trial_timepoint<- factor(MCGILL_df$random_trial_timepoint, levels = c("D0", "D5", "D15", "HLY"))
MCGILL_df$scaled_cfDNA <- 
  ( MCGILL_df$qubit_extracted_cfDNA*MCGILL_df$elution_volume_uL/MCGILL_df$plasma_volume )*
  ( MCGILL_df$batch_nucleic_acid_control_input_concentration*MCGILL_df$batch_nucleic_acid_control_input_volume) /
  ( MCGILL_df$batch_nucleic_acid_control_output * 45)
MCGILL_df$hosp <- MCGILL_df$hospitalization_status_as_of_june_26_2020
MCGILL_df$hosp[grepl("hospita|disch", MCGILL_df$hosp)]<-'Alive'
MCGILL_df$hosp<-factor(MCGILL_df$hosp, levels = c("deceased", "Alive", "HLY"))
MCGILL_df.melt <- melt(MCGILL_df[, c('random_trial_treatment', 'hosp', 'patient_id', 'sample_id', 'scaled_cfDNA', 'random_trial_timepoint', 'sample_class', tissues)], 
                       id.vars = c('random_trial_treatment', 'hosp', 'patient_id', 'sample_id', 'scaled_cfDNA', 'random_trial_timepoint', 'sample_class'))

# Counts -------
# Num deceased
nrow(MCGILL_df[MCGILL_df$hosp=="deceased", ])
nrow(MCGILL_df[MCGILL_df$hosp=="Alive", ])
nrow(MCGILL_df[MCGILL_df$hosp=="HLY", ])

table((MCGILL_df[, c("hosp", "random_trial_timepoint")]))
# Standard of care vs rito/lopi 
pvals <- c()
for (tissue in tissues){
  sub <- MCGILL_df[!grepl("healthy_control", MCGILL_df$cohort), ]
  sub <- melt(sub[, c('random_trial_treatment', 'sample_id', tissues)], id.vars = c("random_trial_treatment", "sample_id"))
  sub <- sub[sub$variable == tissue, ]
  print("================")
  print(tissue)
  print(wilcox.test(sub$value~sub$random_trial_treatment)$p.value)
  pvals <- c(pvals, wilcox.test(sub$value~sub$random_trial_treatment)$p.value)
}

# Dead vs Alive vs Healthy
pvals <- c()
for (tissue in tissues){
  sub <- MCGILL_df.melt[MCGILL_df.melt$variable == tissue, ]
  #sub$value <- sub$value * sub$scaled_cfDNA
  
  print("======================")
  print(tissue)
  print("==PVALUE==")
  print(pairwise.wilcox.test(sub$value, sub$hosp, p.adjust.method = "none"))
  print("==MEAN==")
  print(aggregate(.~hosp, sub[, c('hosp', 'value')], mean))
  #pvals <- c(pvals, wilcox.test(sub2$value~sub2$hosp)$p.value)
}
# COVID vs Healthy
pvals <- c()
for (tissue in tissues){
  sub <- MCGILL_df.melt[MCGILL_df.melt$variable == tissue, ]
  sub$value <- sub$scaled_cfDNA * sub$value
  print("======================")
  print(tissue)
  print(wilcox.test(sub$value~sub$sample_class)$p.value)
  print(aggregate(.~sample_class, sub[, c('sample_class', 'value')], mean))
  pvals <- c(pvals, wilcox.test(sub$value~sub$sample_class)$p.value)
}

# ordinal scale analysis

ordinals <- unique(MCGILL_df[grepl("MCGILL", MCGILL_df$sample_id), grepl("ordinal|patient", colnames(MCGILL_df))])
ordinals <- melt(ordinals, id.vars = c("patient_id"))
ordinals$value <- as.numeric(as.character(ordinals$value))
ordinals$variable<-as.numeric(gsub("ordinal_D", "", ordinals$variable))

ordinals <- ordinals[order(ordinals$value), ]
ggplot(data= ordinals[!is.na(ordinals$value), ]) + geom_point(aes(x=variable, y=value))+
  geom_line(aes(x=variable, y=value))+
  geom_point(data = ordinals2, aes(x=as.numeric(gsub("D", "", random_trial_timepoint)), 
                                   y= scaled_cfDNA*10/max(ordinals2$scaled_cfDNA)), color = 'blue')+
  facet_wrap(vars(patient_id))+
  scale_y_continuous(breaks = 0:10, sec.axis = sec_axis(~./10*max(ordinals2$scaled_cfDNA), name = "Total cfDNA concentration (ng/uL)"))+
  ylab("Ordinal scale score")+
  xlab("Clinical trial timepoint")+
  theme_alex()+
  theme(strip.text = element_text(family = "Helvetica", size = 6),
        axis.title.y.right = element_text(family = "Helvetica", size = 8, color = "blue"))

ordinals2 <- MCGILL_df[!grepl("SWIFT", MCGILL_df$sample_id), ]
ordinals2$ord <- NA
for (i in 1:nrow(ordinals2)){
  if (ordinals2$random_trial_timepoint[i] == "D0"){
    ordinals2$ord[i]<- ordinals2$ordinal_D0[i]
  }
  if (ordinals2$random_trial_timepoint[i] == "D1"){
    ordinals2$ord[i]<- ordinals2$ordinal_D0[i]
  }
  if (ordinals2$random_trial_timepoint[i] == "D5"){
    ordinals2$ord[i]<- ordinals2$ordinal_D5[i]
  }
  if (ordinals2$random_trial_timepoint[i] == "D15"){
    ordinals2$ord[i]<- ordinals2$ordinal_D15[i]
  }
}
ordinals.m <- ordinals2[, c('ord', 'lung', 'erythroblast', 'scaled_cfDNA')]
ordinals.m$lung <- ordinals.m$lung*ordinals.m$scaled_cfDNA
ordinals.m$erythroblast <- ordinals.m$erythroblast*ordinals.m$scaled_cfDNA
#ordinals.m$liver <- ordinals.m$liver*ordinals.m$scaled_cfDNA
ordinals.m <- melt(ordinals.m, id.vars = 'ord')
#ordinals.m$ord <- factor(ordinals.m$ord, levels = c("4", "5", "6", "7", "8", "9"))
ordinals.m$ord <- as.numeric(ordinals.m$ord)
if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/ORDINAL_ALL_CONC.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = ordinals.m)+
  annotate(geom = "rect", xmin = (6.5), xmax = (9.5), ymin = 0, ymax = 6, fill = 'lightblue', alpha = 0.5)+
  geom_boxplot(aes(x=(ord), y=value, color = variable, group = interaction(ord, variable)), outlier.shape = NA, position = position_dodge((width = 0.9)))+
  geom_point(aes(x=(ord), y=value, color = variable),
             position = position_dodge(width = 0.9))+
  xlab("WHO ordinal score")+
  scale_color_manual(values = as.vector(pals::polychrome(19))[c(16,4, 1)])+
  scale_y_continuous(labels = function(x){sprintf("%.1f", x)}, limits = c(-0.25, 6))+
  scale_x_continuous(breaks = c(4:9))+
  ylab("cfDNA concentration (ng/uL)")+
  theme_alex()+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig)dev.off()

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/ORDINAL_ALL_CONC_ZOOM.pdf",
                 width=29.5/25.4, height=25/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = ordinals.m[grepl("lung|eryth", ordinals.m$variable) & grepl("4|5|6|7|8", ordinals.m$ord), ])+
  annotate(geom = "rect", xmin = (6.5), xmax = (8.5), ymin = 0, ymax = 0.35, fill = 'lightblue', alpha = 0.5)+
  geom_boxplot(aes(x=(ord), y=value, color = variable, group = interaction(ord, variable)), outlier.shape = NA, position = position_dodge((width = 0.9)))+
  geom_point(aes(x=(ord), y=value, color = variable),
             position = position_dodge(width = 0.9),
             size = 1, pch=21)+
  xlab("Clinical progression score")+
  scale_color_manual(values = as.vector(pals::polychrome(19))[c(16,4, 1)])+
  ylab("cfDNA concentration (ng/uL)")+
  theme_alex()+
  theme(plot.margin = unit(c(2,2,2,2), "pt"),
        axis.title = element_blank())
if(save_fig)dev.off()

ordinals.m$group <- NA
ordinals.m$group[as.numeric(as.character(ordinals.m$ord)) >= 7]<-'7-9 (n=22)'
ordinals.m$group[as.numeric(as.character(ordinals.m$ord)) <7 ]<-'4-6 (n=30)'


oo <- ordinals.m[ordinals.m$variable == "erythroblast", ]
wilcox.test(oo$value~oo$group)
t.test(oo$value~oo$group)
a<-pROC::roc(oo$group, oo$value, ci=TRUE, plot=TRUE, print.auc=TRUE)
cor.test(oo$ord, oo$value, method = 'spearman')

oo2 <- oo[grepl("7|8|9", oo$ord), ]
oo2$ord<- as.character(oo2$ord)
oo2$ord[oo2$ord=="8"]<-"7"
t.test(oo2$value~oo2$ord)
wilcox.test(oo2$value~oo2$ord)

oo <- ordinals.m[ordinals.m$variable == "lung", ]
wilcox.test(oo$value~oo$group)
t.test(oo$value~oo$group)
b<-pROC::roc(oo$group, oo$value, ci=TRUE, plot=TRUE, print.auc=TRUE)
cor.test(oo$ord, oo$value, method = 'pearson')
ggplot(data = oo)+
  geom_point(aes(x=ord, y = value))
oo$over <- FALSE
oo$over[oo$value > 0.005]<-TRUE
aggregate(.~group, oo[, c('group', 'over')], sum)
table(oo$group)

prop.test(c(3, 10), c(30, 22))

oo <- ordinals.m[ordinals.m$variable == "scaled_cfDNA", ]
t.test(oo$value~oo$group)
wilcox.test(oo$value~oo$group)
c<-pROC::roc(oo$group, oo$value, ci=TRUE, plot=TRUE, print.auc=TRUE)
cor.test(oo$ord, oo$value, method = 'spearman')


roc <- data.frame(c(a$sensitivities, b$sensitivities, c$sensitivities), 
                  c(a$specificities, b$specificities, c$specificities), 
                  c(rep("ery", length(a$sensitivities)), rep("lung", length(b$sensitivities)), rep("total", length(c$sensitivities))))



colnames(roc)<-c('sensitivity', 'specificity', 'variable')

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/ORDINAL_ROC_CONC.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = roc)+
  geom_abline(slope = 1, intercept = 0, color = "gray")+
  geom_path(aes(x=1-specificity, y= sensitivity, color = variable))+
  scale_color_manual(values = as.vector(pals::polychrome(19))[c(4,16, 1)])+
  theme_alex()+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig)dev.off()


# FIGURES ------------------------------------------------------------------------------------
if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_frac_byday_mega.pdf",
                 width=89/25.4, height=89/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = MCGILL_df.melt[grepl("ery|neut|lung|liver", MCGILL_df.melt$variable), ])+
  geom_boxplot(aes(x=random_trial_timepoint, y= value, color = variable, fill = hosp), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=random_trial_timepoint, y= value, group = interaction(hosp, variable), color = variable, shape = hosp),
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c(NA, NA, NA))+
  scale_color_manual(values = as.vector(pals::polychrome(19))[c(9,4,16,17)])+
  facet_wrap(vars(variable), scales = "free_x", ncol = 2)+
  theme_alex()+ylab("Tissue fraction")+xlab("Trial timepoint")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_conc_byday_mega.pdf",
                 width=89/25.4, height=89/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = MCGILL_df.melt[grepl("ery|neut|lung|liver", MCGILL_df.melt$variable), ])+
  geom_boxplot(aes(x=random_trial_timepoint, y= value*scaled_cfDNA, color = variable, fill = hosp), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=random_trial_timepoint, y= value*scaled_cfDNA, group = interaction(hosp, variable), color = variable, shape = hosp),
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c(NA, NA, NA))+
  scale_color_manual(values = as.vector(pals::polychrome(19))[c(9,4,16,17)])+
  facet_wrap(vars(variable), scales = "free_y", ncol = 2)+
  theme_alex()+ylab("Tissue fraction")+xlab("Trial timepoint")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

ggplot(data = MCGILL_df.melt[grepl("ery|neut|lung|liver", MCGILL_df.melt$variable), ])+
  geom_density(aes(x=value, color = hosp))+
  facet_wrap(vars(variable), scales = "free")


if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_frac_death_mega.pdf",
                 width=60/25.4, height=50/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = MCGILL_df.melt[grepl("ery|neut|lung|liver", MCGILL_df.melt$variable), ])+
  geom_boxplot(aes(x=variable, y= value, color = hosp), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=variable, y= value, color = hosp),
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c(NA, NA, NA))+
  scale_color_manual(values = c('red', 'purple', 'forestgreen'))+
  scale_y_continuous(labels = function(x){sprintf("%.2f", x)})+
  theme_alex()+ylab("Tissue fraction")+xlab("Tissue")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_frac_death_all_tissues_mega.pdf",
                 width=240/25.4, height=70/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = MCGILL_df.melt)+
  geom_boxplot(aes(x=variable, y= value, color = hosp), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=variable, y= value, color = hosp),
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c(NA, NA, NA))+
  scale_color_manual(values = c('red', 'purple', 'forestgreen'))+
  theme_alex()+ylab("Tissue fraction")+xlab("Tissue")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_frac_std_vs_exp_mega.pdf",
                 width=240/25.4, height=192/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = MCGILL_df.melt[MCGILL_df.melt$random_trial_treatment!="HLY", ])+
  geom_boxplot(aes(x=random_trial_treatment, y= value), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=random_trial_treatment, y= value),
             position = position_dodge(width = 0.9))+
  theme_alex()+ylab("Tissue fraction")+xlab("Tissue")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  facet_wrap(vars(variable))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_livervsAST_mega.pdf",
                 width=30/25.4, height=30/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=MCGILL_df[!grepl("HLY", MCGILL_df$AST), ])+
  annotate(geom = "rect", xmin = 10, xmax = 40, ymin = 0, ymax = 0.45, fill = 'forestgreen', alpha = 0.5)+
  #geom_smooth(aes(x=as.numeric(as.character(AST)), y = liver), method = 'lm', se = FALSE)+
  geom_point(aes(x=as.numeric(as.character(AST)), y = liver))+
  theme_alex()+
  xlab('Aspartate aminotransferase (units/L)')+
  ylab('Liver fraction')+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}
if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_livervsALT_mega.pdf",
                 width=30/25.4, height=30/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=MCGILL_df)+
  #geom_smooth(aes(x=as.numeric(as.character(AST)), y = liver), method = 'lm', se = FALSE)+
  annotate(geom = "rect", xmin = 7, xmax = 56, ymin = 0, ymax = 0.45, fill = 'forestgreen', alpha = 0.5)+
  geom_point(aes(x=as.numeric(as.character(ALT)), y = liver))+
  theme_alex()+
  xlab('(units/L)')+
  ylab('Liver fraction')+
  scale_x_log10()+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}
if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_eryth_vs_hematocrit_mega.pdf",
                 width=30/25.4, height=30/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=MCGILL_df)+
  annotate(geom = "rect", xmin = 0.36, xmax = 0.5, ymin = 0, ymax = 0.7, fill = 'forestgreen', alpha = 0.5)+
  #geom_smooth(aes(x=as.numeric(as.character(hematocrit)), y = erythroblast), method = 'lm', se = FALSE)+
  geom_point(aes(x=as.numeric(as.character(hematocrit)), y = erythroblast))+
  theme_alex()+
  xlab('Hematocrit')+
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5))+
  ylab('Erythroblast fraction')+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_eryth_vs_hemo_mega.pdf",
                 width=30/25.4, height=30/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=MCGILL_df)+
  #geom_smooth(aes(x=as.numeric(as.character(hematocrit)), y = erythroblast), method = 'lm', se = FALSE)+
  annotate(geom = "rect", xmin = 12, xmax = 17.5, ymin = 0, ymax = 0.7, fill = 'forestgreen', alpha = 0.5)+
  geom_point(aes(x=as.numeric(as.character(hemoglobin)), y = erythroblast))+
  theme_alex()+
  xlab('Hemoglobin')+
  ylab('Erythroblast fraction')+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

aa <- MCGILL_df[MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA>0.1, ]

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/MCGILL_total cfDNA_vs_LDH_mega.pdf",
                 width=30/25.4, height=30/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6)}
ggplot(data=MCGILL_df)+
  #geom_smooth(aes(x=as.numeric(as.character(LDH)), y = scaled_cfDNA), method = 'lm', se = FALSE)+
  annotate(geom = "rect", xmin = 140, xmax = 280, ymin = 0, ymax = 6, fill = 'forestgreen', alpha = 0.5)+
  geom_point(aes(x=as.numeric(as.character(LDH)), y = scaled_cfDNA))+
  theme_alex()+
  xlab('Lactate dehydrogenase (units/L)')+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x))+
  ylab('Total cfDNA concentration')+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

ggplot(data=MCGILL_df)+
  #geom_smooth(aes(x=as.numeric(as.character(LDH)), y = scaled_cfDNA), method = 'lm', se = FALSE)+
  #annotate(geom = "rect", xmin = 140, xmax = 280, ymin = 0, ymax = 6, fill = 'forestgreen', alpha = 0.5)+
  geom_point(aes(x=as.numeric(as.character(LDH)), y = erythroblast*scaled_cfDNA))+
  theme_alex()+
  xlab('Lactate dehydrogenase (units/L)')+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x))+
  ylab('Total cfDNA concentration')+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))

cor.test(as.numeric(MCGILL_df$hematocrit), MCGILL_df$erythroblast, method = 'pearson')
cor.test(as.numeric(MCGILL_df$hematocrit), MCGILL_df$erythroblast, method = 'spearman')

cor.test(as.numeric(MCGILL_df$hemoglobin), MCGILL_df$erythroblast, method = 'pearson')
cor.test(as.numeric(MCGILL_df$hemoglobin), MCGILL_df$erythroblast, method = 'spearman')

cor.test(as.numeric(MCGILL_df$ALT), MCGILL_df$liver, method = 'pearson')
cor.test(as.numeric(MCGILL_df$ALT), MCGILL_df$liver, method = 'spearman')

cor.test(as.numeric(MCGILL_df$AST), MCGILL_df$liver, method = 'pearson')
cor.test(as.numeric(MCGILL_df$AST), MCGILL_df$liver, method = 'spearman')

cor.test(as.numeric(MCGILL_df$creatinine), MCGILL_df$kidney, method = 'pearson')
cor.test(as.numeric(MCGILL_df$creatinine), MCGILL_df$kidney, method = 'spearman')

cor.test(as.numeric(MCGILL_df$hematocrit), MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA, method = 'pearson')
cor.test(as.numeric(MCGILL_df$hematocrit), MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA, method = 'spearman')

cor.test(as.numeric(MCGILL_df$hemoglobin), MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA, method = 'pearson')
cor.test(as.numeric(MCGILL_df$hemoglobin), MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA, method = 'spearman')

cor.test(as.numeric(MCGILL_df$ALT), MCGILL_df$liver*MCGILL_df$scaled_cfDNA, method = 'pearson')
cor.test(as.numeric(MCGILL_df$ALT), MCGILL_df$liver*MCGILL_df$scaled_cfDNA, method = 'spearman')

cor.test(as.numeric(MCGILL_df$AST), MCGILL_df$liver*MCGILL_df$scaled_cfDNA, method = 'pearson')
cor.test(as.numeric(MCGILL_df$AST), MCGILL_df$liver*MCGILL_df$scaled_cfDNA, method = 'spearman')

cor.test(as.numeric(MCGILL_df$creatinine), MCGILL_df$kidney*MCGILL_df$scaled_cfDNA, method = 'pearson')
cor.test(as.numeric(MCGILL_df$creatinine), MCGILL_df$kidney*MCGILL_df$scaled_cfDNA, method = 'spearman')

cor.test(as.numeric(MCGILL_df$LDH), MCGILL_df$erythroblast, method = "pearson")
cor.test(as.numeric(MCGILL_df$LDH), MCGILL_df$erythroblast, method = "spearman")

cor.test(as.numeric(MCGILL_df$LDH), MCGILL_df$scaled_cfDNA, method = "pearson")
cor.test(as.numeric(MCGILL_df$LDH), MCGILL_df$scaled_cfDNA, method = "spearman")



sub <- MCGILL_df[MCGILL_df$sample_class=="SARS-CoV2+", ]
pROC::roc(MCGILL_df$hosp, MCGILL_df$erythroblast, ci= TRUE)
pROC::roc(MCGILL_df$hosp, MCGILL_df$erythroblast*MCGILL_df$scaled_cfDNA, ci = TRUE)
pROC::roc(MCGILL_df$hosp, MCGILL_df$scaled_cfDNA, ci = TRUE)

# TESTING -----

AA <- MCGILL_df

AA$lung_tot <- AA$lung*AA$scaled_cfDNA

AA <- AA[, c('patient_id', 'sample_class', 'lung_tot', 'hosp')]

AA$lotsoflung <- AA$lung_tot>8.149*10**-4
aggregate(.~sample_class, AA[, c('lotsoflung', 'sample_class')], sum)
table(AA$sample_class)

prop.test(c(0,32), c(4,52), alternative = 'two.sided', correct = FALSE)

AA <- MCGILL_df

AA <- AA[, c('patient_id', 'sample_class', 'erythroblast', 'hosp')]

AA$loteryt <- AA$erythroblast>0.0
aggregate(.~sample_class, AA[, c('loteryt', 'sample_class')], sum)
table(AA$sample_class)

prop.test(c(0,34), c(4,52), alternative = 'two.sided', correct = FALSE)
