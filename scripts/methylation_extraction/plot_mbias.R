library(data.table)
library(ggplot2)
source('~/theme_alex_minimal.R')
args = commandArgs(trailingOnly=TRUE)

seq_type = args[[1]]
out=args[[2]]
df = fread('cat /dev/stdin', fill=TRUE, na.strings = 'NA')
colnames(df) <-c('pos', 'count_meth', 'count_unmeth', 'pmeth', 'cov')
df$pmeth <- as.numeric(as.character(df$pmeth))

df$context <- NA
df$read <-NA
df$trim_status <- NA

cytosines <- c('CpG', 'CHG', 'CHH')

c_count <- 0
r_count <- 0
t_count <- 0

if (grepl('1x', seq_type)){
  sets <- rep(cytosines, 2) #trim and untrimmed
  reads <- c('R1', 'R1') #trim and untrimmed
  trim <- c('trim', 'untrimmed')
  
  for (i in 1:nrow(df)){
    if (df$pos[i] == 1){
      c_count <- c_count + 1
    }
    if (c_count == 4){
      t_count <- t_count + 1
    }
    df$context <- sets[c_count]
    df$read[i] <- 'R1' #because single-end
    df$trim_stats <- trim[t_count]
  }
  
}

if (grepl('2x', seq_type)){
  sets <- rep(cytosines, 4) #R1, R2, trimmed & untrimmed
  reads <- c('R1', 'R2', 'R1', 'R2')
  trim <- c('trim', 'untrimmed')
  
  for (i in 1:nrow(df)){
    if (df$pos[i] == 1 ){
      c_count <- c_count + 1
    }
    if (c_count == 4 | c_count == 7 | c_count == 10){
      r_count <- r_count + 1
    }
    if (r_count == 3){
      t_count <- t_count + 1
    }
    df$context <- sets[c_count]
    df$read[i] <- reads[r_count] #because single-end
    df$trim_stats <- trim[t_count]
  }
}

df$context<-factor(df$context, levels=c('CpG', 'CHG', 'CHH'))

pdf(file =out, width = 6, height=6, pointsize = 8)
ggplot(data=df)+
  geom_line(aes(x=pos, y= pmeth, color=read))+
  facet_grid(context~trim_status, scales='free_y')+
  xlab('position')+
  ylab('percent methylation')+
  scale_color_manual(values=c('steelblue', 'orange'))+
  theme(legend.position = "top",
        legend.title = element_text(family="Helvetica", size=8),
        legend.text = element_text(family="Helvetica", size=8))+
  theme_alex_minimal()
  

dev.off()
