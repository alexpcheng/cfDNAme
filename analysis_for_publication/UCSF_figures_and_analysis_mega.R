#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Initialize workspace -----------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(ggplot2)
source('~/theme_alex.R')
master_table_path <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/master_table_megakaryocytes.csv'
figure_path <- '/workdir/apc88/WGBS_pipeline/Figures/'
workdir <- '/workdir/apc88/WGBS_pipeline/analysis_for_publication/'
setwd(workdir)

# Functions -----------------------------------------------------------------------------------
hour_to_decimal <- function(timestr){
  hour <- as.numeric(strsplit(x=timestr, split = ':')[[1]][1])
  minutes <- as.numeric(strsplit(x=timestr, split = ':')[[1]][2])
  hour <- (hour + minutes/60)/24
  return(hour)
}

# Read in master table -----------------------------------------------------------------------
df <- fread(master_table_path, na.strings = "NA", data.table = FALSE)

tissues <- c('macrophage', 'monocyte', 'dendritic','neutrophil', 'eosinophil', 'Tcell', 'Nkcell', 'Bcell',
             'erythroblast', 'hsc', 'progenitor','spleen', 'kidney', 'heart', 'skin', 'lung','liver','pancreas',
             'colon', 'megakaryocyte')
# Longitudinal analysis of COVID19 -----------------------------------------------------------

UCSF_df <- df[df$cohort == "high_frequency_collection", ]
UCSF_df$sample_class <- factor(UCSF_df$sample_class, levels = c("SARS-CoV2+", "other_RNA_virus"))
UCSF_df$date <- as.Date(gsub(" .*", "", UCSF_df$`days_past_covid_symptom_onset_normalized_01/01/00`), format = "%m/%d/%y")
UCSF_df$time <- sapply(X = gsub(".* ", "", UCSF_df$`days_past_covid_symptom_onset_normalized_01/01/00`),
                       FUN = hour_to_decimal)

UCSF_df$day <- as.numeric((UCSF_df$date - as.Date("01/01/00", format = "%m/%d/%y")+UCSF_df$time))
UCSF_df$scaled_cfDNA <- UCSF_df$hg19_aligned_reads/UCSF_df$spikein_aligned_reads/UCSF_df$plasma_volume

UCSF_df.melt <- melt(UCSF_df[, c('patient_id', 'sample_id', 'sample_class', 'cohort', tissues, 'day', 'scaled_cfDNA')], id.vars = c('patient_id', 'sample_id', 'sample_class',
                                                                'cohort', 'day', 'scaled_cfDNA'))

UCSF_df.melt$variable<- factor(UCSF_df.melt$variable, levels=tissues)


for (tissue in tissues){
  sub <- UCSF_df.melt[UCSF_df.melt$variable == tissue, c('sample_class', 'value')]
  print("=======================================")
  print(tissue)
  print(aggregate(.~sample_class, sub[, c('value', 'sample_class')], mean))
  print(wilcox.test(sub$value~sub$sample_class)$p.value)
}

# TIME BETWEEN COLLECTIONS -------------------------------------------------------------------
time_df <- UCSF_df[UCSF_df$sample_class == "SARS-CoV2+", ]
table(time_df$patient_id)

time_df <- time_df[order(time_df$patient_id, time_df$day), c('patient_id', 'day')]

time_df <- time_df %>% group_by(patient_id) %>% mutate(newdiff = day - lag(day))
mean(time_df$newdiff, na.rm = TRUE)

# DISSIMILARITY ------------------------------------------------------------------------------
# Dissimilarity between patients, days, and timepoints
dissimi <- UCSF_df[UCSF_df$sample_class=="SARS-CoV2+", c('patient_id', 'sample_id', tissues, 'day')]
bray_curtis <- data.frame(as.matrix((vegan::vegdist(x=dissimi[, tissues], method = "bray"))))
colnames(bray_curtis) <- dissimi$sample_id
rownames(bray_curtis) <- dissimi$sample_id

a <- bray_curtis

b <- melt(as.matrix(a))
b <- b[!is.na(b$value), ]
colnames(b)<- c('sample_id', 'sample_id2', 'value')
b <- merge(b, dissimi[, c('sample_id', 'patient_id', 'day')]) 
colnames(b)<-c('sample_id1', 'sample_id', 'value', 'patient_id1', 'day1')
b <- merge(b, dissimi[, c('sample_id', 'patient_id', 'day')])
colnames(b)<-c('sample_id2', 'sample_id1', 'value', 'patient_id1', 'day1', 'patient_id2', 'day2')

b$diff <- abs(b$day1-b$day2)
c <- b[order(b$diff, b$value), ]
c <- c[c$patient_id1 == c$patient_id2, ]
c <- c[, c('diff', 'value', 'patient_id1', 'patient_id2', 'day1', 'day2', 'sample_id1', 'sample_id2')]
c$diff2 <- abs( (c$day1 - floor(c$day1)) - (c$day2 - floor(c$day2)) )
q <- c

q$un <- apply(X=q, MARGIN=1, FUN=
                function(x){
                  a <- x['sample_id1']
                  b <- x['sample_id2']
                  print(a)
                  c<- sort(c(a,b))
                  d <- paste0(c[1], c[2])
                  print(d)
                  return(d)
                })
q <- q[order(q$day1, decreasing=TRUE), ]
q <- q[!duplicated(q$un), ]

r <- b
r$un <- apply(X=r, MARGIN=1, FUN=
                function(x){
                  a <- x['sample_id1']
                  b <- x['sample_id2']
                  print(a)
                  c<- sort(c(a,b))
                  d <- paste0(c[1], c[2])
                  print(d)
                  return(d)
                })
r <- r[!duplicated(r$un), ]
r <- r[r$sample_id1!=r$sample_id2, ]

all_p<-r[r$patient_id1!=r$patient_id2, ]$value
#same_p <- r[r$patient_id1==r$patient_id2 & abs((r$day1) - (r$day2))>=1, ]$value
same_p <- r[r$patient_id1==r$patient_id2, ]$value
same_d <- r[r$patient_id1==r$patient_id2 & abs((r$day1) - (r$day2))<1, ]$value

wilcox.test(same_p, all_p)
wilcox.test(same_p, all_p)
wilcox.test(same_p, all_p)


AA <- UCSF_df
AA <- AA[, c('patient_id', 'sample_class', 'erythroblast')]

AA$lotsofer <- AA$erythroblast>0.00
aggregate(.~sample_class, AA[, c('lotsofer', 'sample_class')], sum)
table(AA$sample_class)

prop.test(c(1,39), c(6,52), alternative = 'two.sided', correct = FALSE)

ggplot(data = UCSF_df)+
  geom_density(aes(x=liver))

sum(UCSF_df$neutrophil[UCSF_df$sample_class=="other_RNA_virus"]<0.01)/6*100
sum(UCSF_df$neutrophil[UCSF_df$sample_class=="other_RNA_virus"]<0.01)/6*100


### ZERO INFLATION ----
sum(UCSF_df$liver[UCSF_df$sample_class=="SARS-CoV2+"]<=0.0)/52*100
sum(UCSF_df$liver[UCSF_df$sample_class=="other_RNA_virus"]<=0.0)/6*100
sum(UCSF_df$liver<=0.0)/58*100

save_fig = T
# FIGURES ------------------------------------------------------------------------------------
if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UCSF_abundance_pass_frac_mega.pdf",
                width=90/25.4, height=60/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot(data=UCSF_df.melt[UCSF_df.melt$sample_class=="SARS-CoV2+", ])+
  geom_area(aes(x=day, y= value, fill=variable))+
  geom_point(data = unique(UCSF_df.melt[UCSF_df.melt$sample_class=="SARS-CoV2+", c('day', 'patient_id')]),
             aes(x=day, y=0), pch=17)+
  scale_fill_manual(values=as.vector(pals::polychrome())[c(1:3,9,5:8,4,10,20,12:19,28)])+
  facet_wrap(vars(patient_id), scales='free_x', ncol=1, strip.position = "right")+
  theme_alex()+
  scale_y_continuous(breaks = c(0,1))+
  ylab('Tissue fraction')+
  theme(panel.spacing = unit(0, 'pt'),
        strip.text = element_text(family='Helvetica', size=6),
        plot.margin=grid::unit(c(1,1,1,1), "pt"))
if(save_fig)dev.off()

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UCSF_abundance_pass_FLUCTL_frac_mega.pdf",
                width=25/25.4, height=46/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot(data=UCSF_df.melt[UCSF_df.melt$sample_class=="other_RNA_virus", ])+
  geom_col(aes(x=patient_id, y= value, fill=variable), color='black')+
  scale_fill_manual(values=as.vector(pals::polychrome(28))[c(1:3,9,5:8,4,10,20,12:19,28)])+
  theme_alex()+
  ylab('tissue fraction')+
  xlab('Influenza positive')+
  theme(axis.title.y = element_blank(),
        plot.margin=grid::unit(c(1,1,1,1), "pt"))
if(save_fig)dev.off()

ggplot(data=UCSF_df.melt[grepl("eryth|liver|neutrophil|lung", UCSF_df.melt$variable), ])+
  geom_density(aes(x=value))+
  facet_wrap(vars(variable), scales = "free")

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UCSF_frac_mega.pdf",
                width=60/25.4, height=35/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot(data=UCSF_df.melt[grepl("eryth|liver|neutrophil|lung", UCSF_df.melt$variable), ])+
  geom_boxplot(aes(x=variable, y= value, color=sample_class), pch=21, fill=NA, outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=variable, y= value, color = sample_class), position = position_dodge(width=0.9))+
  scale_color_manual(values = c('steelblue', 'orange', 'forestgreen'))+
  theme_alex()+#+ylim(c(0,0.3))+
  xlab('Tissue')+
  ylab('Tissue fraction')+
  theme(plot.margin=grid::unit(c(2,1,1,1), "pt"))
if(save_fig)dev.off()

if(save_fig){pdf(file="/workdir/apc88/WGBS_pipeline/Figures/UCSF_frac_death_all_tissues_mega.pdf",
                 width=240/25.4, height=70/25.4, paper="special", bg="white",
                 fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats = FALSE)}
ggplot(data = UCSF_df.melt)+
  geom_boxplot(aes(x=variable, y= value, color = sample_class), outlier.shape = NA, width = 0.9)+
  geom_point(aes(x=variable, y= value, color = sample_class),
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values= c(NA, NA, NA))+
  scale_color_manual(values = c('steelblue', 'orange'))+
  theme_alex()+ylab("Tissue fraction")+xlab("Tissue")+
  theme(strip.text = element_text(family = "Helvetica", size = 6,
                                  margin = margin(2,0,2,0,"pt")))+
  theme(plot.margin = unit(c(2,2,2,2), "pt"))
if(save_fig){dev.off()}

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/Bray_by_day_mega.pdf",
                width=45/25.4, height=35/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot(data=q[!duplicated(q$un) & q$sample_id1!=q$sample_id2, ])+
  geom_point(aes(x=diff, y=value, color = patient_id1), pch=21)+
  geom_smooth(aes(x=diff, y=value, color = patient_id1), method = 'lm', se = FALSE)+
  xlab("Time between two collections (days)")+
  ylab("Bray-Curtis dissimilarity")+
  scale_y_continuous(limits=c(0,0.73))+
  ggsci::scale_color_nejm()+
  theme_alex()+
  theme(plot.margin=grid::unit(c(2,2,1,1), "pt"))
if(save_fig)dev.off()

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/Bray_by_class_mega.pdf",
                width=16/25.4, height=35/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot()+
  geom_boxplot(data = r[r$patient_id1!=r$patient_id2, ], aes(x='cAll samples', y=value))+
  geom_boxplot(data = r[r$patient_id1==r$patient_id2 & abs((r$day1) - (r$day2))>=1, ], aes(x='bSame patient', y=value))+
  geom_boxplot(data= r[r$patient_id1==r$patient_id2 & abs((r$day1) - (r$day2))<1, ], aes(x='aSame day', y=value))+
  geom_point(data = r[r$patient_id1!=r$patient_id2, ], aes(x='cAll samples', y=value))+
  geom_point(data = r[r$patient_id1==r$patient_id2 & floor(r$day1) != floor(r$day2), ], aes(x='bSame patient', y=value))+
  geom_point(data= r[r$patient_id1==r$patient_id2 & floor(r$day1) == floor(r$day2), ], aes(x='aSame day', y=value))+
  theme_alex()+scale_y_continuous(limits=c(0,0.73))+
  ylab("Bray Curtis dissimilarity")+
  xlab(" ")+
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin=grid::unit(c(2,1,1,1), "pt"))
if(save_fig)dev.off()

if(save_fig)pdf(file="/workdir/apc88/WGBS_pipeline/Figures/Bray_by_heatmap_mega.pdf",
                width=165/25.4, height=30/25.4, paper="special", bg="white",
                fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)
ggplot(data = q)+
  geom_tile(aes(x=factor(round(day1, digits = 1)), y= factor(round(day2, digits = 1)), fill = (value)), color = 'white')+
  facet_wrap(vars(patient_id1), scales = "free", nrow=1)+
  xlab("day")+
  ylab("day")+
  theme_alex()+
  theme(panel.grid = element_blank(), strip.text = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        plot.margin=grid::unit(c(1,1,1,1), "pt"),
        axis.title = element_blank())+
  scale_fill_viridis_c()
if(save_fig)dev.off()

