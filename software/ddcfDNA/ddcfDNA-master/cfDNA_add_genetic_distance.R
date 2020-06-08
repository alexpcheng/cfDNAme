######################################################################
# getting input parameters
######################################################################

# if you want see commands in output file
# options(echo=TRUE) 

# getting the arguments
args <- commandArgs(trailingOnly = TRUE)

# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))
print(args)

######################################################################
# loading librares
######################################################################

if (!require("reshape2")) {
  library(reshape2)
}


#####################################################################-
# util functions
######################################################################

#####################################################################-
# Setting input parameters
######################################################################

# setting the value of the arguments
genotype_filename <- args[1]
decode_dir <- args[2]
output_filename <- args[3]


print(sprintf("genotype_filename:%s",genotype_filename))
print(sprintf("decode_dir:%s",decode_dir))
print(sprintf("output_filename:%s",output_filename))


######################################################################-
# loading tables
#######################################################################

geno_df = read.delim(genotype_filename)

gmap_male_nc = read.delim( paste(decode_dir, "/", "male_noncarrier.gmap.noSpaces" ,sep="") )
gmap_female_nc = read.delim( paste(decode_dir, "/", "female_noncarrier.gmap.noSpaces" ,sep="") )

gmap_male_c = read.delim( paste(decode_dir, "/", "male_carrier.gmap.noSpaces" ,sep="") )
gmap_female_c = read.delim( paste(decode_dir, "/", "female_carrier.gmap.noSpaces" ,sep="") )

rmap_male = read.delim( paste(decode_dir, "/", "male.rmap.noSpaces" ,sep="") )
rmap_female = read.delim( paste(decode_dir, "/", "female.rmap.noSpaces" ,sep="") )
rmap_nosex = read.delim( paste(decode_dir, "/", "sex-averaged.rmap.noSpaces" ,sep="") )


######################################################################-
# compute genetic disctance for snps in the data
#######################################################################


# compute comulative distance in decode

gmap_male_nc$cM_num = as.numeric(as.character(gmap_male_nc$cM))
gmap_female_nc$cM_num = as.numeric(as.character(gmap_female_nc$cM))

m_dist_df = gmap_male_nc[is.finite(gmap_male_nc$cM_num),]
f_dist_df = gmap_female_nc[is.finite(gmap_female_nc$cM_num),]

#gmap_male_nc$cM_num_cumsum = cumsum(gmap_male_nc$cM_num)
#gmap_female_nc$cM_num_cumsum = cumsum(gmap_female_nc$cM_num)

rmap_nosex$stdrate_num = as.numeric(as.character(rmap_nosex$stdrate))
ns_rec_df =  rmap_nosex[is.finite(rmap_nosex$stdrate_num ),]



# sort genotyping data
geno_df$ChrNumeric = as.numeric(as.character(substring(geno_df$Chr, first=4)))
geno_df = geno_df[with(geno_df,order(ChrNumeric,Position)),]

geno_df$genetic_dist_log = 0.0
geno_df$genetic_dist_cumsum = 0.0
geno_df$genetic_rec_hotspot = 0.0
geno_df$genetic_block_breaks = 0.0
geno_df$genetic_blocks = 0.0
geno_df$genetic_avg_distance_FromPrevBlock = 0.0

chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
  "chr21","chr22")

block_cM_size = 2
cur_max_black_num = 0

for (ch in 1:length(chrs)) {
  
  cur_chr = chrs[ch]
  
  cur_geno_inds = which(geno_df$Chr == cur_chr)
  
  # extereme case in which a chromosome has zero blocks
  if (length(cur_geno_inds)<1) {
    next
  }
  
  cur_m_dist_inds = which(m_dist_df$chr == cur_chr)
  cur_f_dist_inds = which(f_dist_df$chr == cur_chr)
  
  # compute genetic distance
  
  cur_ch_m_dist = m_dist_df[cur_m_dist_inds,]
  cur_ch_f_dist = f_dist_df[cur_f_dist_inds,]
  
  cur_ch_m_dist$cM_num_cumsum = cumsum(cur_ch_m_dist$cM_num)
  cur_ch_f_dist$cM_num_cumsum = cumsum(cur_ch_f_dist$cM_num)
  
  cur_geno_pos = as.numeric(geno_df$Position[cur_geno_inds])
  
  cur_m_geno_cumsum =   approx(cur_ch_m_dist$pos,cur_ch_m_dist$cM_num_cumsum,cur_geno_pos,method = "linear")
  cur_f_geno_cumsum =   approx(cur_ch_f_dist$pos,cur_ch_f_dist$cM_num_cumsum,cur_geno_pos,method = "linear")
 
  cur_m_geno_dist = c(Inf, diff(cur_m_geno_cumsum$y))
  cur_f_geno_dist = c(Inf, diff(cur_f_geno_cumsum$y))
  
  cur_m_geno_dist[is.na(cur_m_geno_dist)] = 0
  cur_f_geno_dist[is.na(cur_f_geno_dist)] = 0
  
  
  cur_avg_geno_dist = (cur_m_geno_dist+cur_f_geno_dist)/2
  
  # extreme case in which there is only one block in a chromosome
  if (length(cur_avg_geno_dist)<2) {
    cur_avg_geno_dist_cumsum = vector(mode="numeric", length=length(cur_avg_geno_dist))
    cur_avg_geno_dist_cumsum[] = Inf
    
    block_avg_dist_cumsum = vector(mode="numeric", length=length(cur_avg_geno_dist))
    block_avg_dist_cumsum[] = Inf
    
    avg_distance_FromPrevBlock = vector(mode="numeric", length=length(cur_avg_geno_dist))
    avg_distance_FromPrevBlock[] = Inf

  } else {
    cur_avg_geno_dist_cumsum = c(0,cumsum(cur_avg_geno_dist[2:length(cur_avg_geno_dist)]))
  
  
  
    # compute recombinations hotspots
    
    cur_ns_rec_inds = which(ns_rec_df$chr == cur_chr)
    cur_ch_ns_rec = ns_rec_df[cur_ns_rec_inds,]
    cur_ch_ns_rec$hotspot = 1.0 * as.numeric(cur_ch_ns_rec$stdrate_num>=10)
    cur_ch_ns_rec$hotspot[is.na(cur_ch_ns_rec$hotspot)] = 0
    
    cur_ns_rec_hotspot_cumsum =   approx(cur_ch_ns_rec$pos,cumsum(cur_ch_ns_rec$hotspot),cur_geno_pos,method = "linear")
    #cur_ns_rec_hotspot_cumsum$y[is.na(cur_ns_rec_hotspot_cumsum$y)] = 0
    
    cur_ns_rec_hotspot_num = c(0,diff(cur_ns_rec_hotspot_cumsum$y))
    
    cur_ns_rec_hotspot_num[is.na(cur_ns_rec_hotspot_num)] = 0
  
    # spliting into blocks
    cur_avg_geno_blocks = vector(mode= "numeric", length(cur_avg_geno_dist))
    cur_avg_geno_blocks_breaks = vector(mode= "numeric", length(cur_avg_geno_dist))
    
    prev_hotspot_cnt = Inf
    
    for (s in 1:(length(cur_avg_geno_dist_cumsum)-2))
    {
      cur_dist_cs = cur_avg_geno_dist_cumsum[s]
      cur_dist_cs_next = cur_avg_geno_dist_cumsum[s+1]
      
      cur_hotspot_cnt = cur_ns_rec_hotspot_cumsum$y[s]
      
      if (is.infinite(prev_hotspot_cnt)) {
        prev_hotspot_cnt = cur_ns_rec_hotspot_cumsum$y[s]
      }
      
      if (s==1 || 
         (cur_dist_cs-prev_dist_cs) >=block_cM_size  || 
         ( (cur_dist_cs_next-prev_dist_cs) >= block_cM_size &&  (cur_dist_cs_next-prev_dist_cs)-block_cM_size -  (block_cM_size-(cur_dist_cs-prev_dist_cs)) > 0  ) )
        #(is.finite(cur_hotspot_cnt-prev_hotspot_cnt) && cur_hotspot_cnt-prev_hotspot_cnt >0.99) )
      {
        if (s>1){
          print(sprintf("%s SELECTED %d, %f --- %f, %f, %f",cur_chr,s, cur_dist_cs-prev_dist_cs,cur_dist_cs_next-prev_dist_cs, cur_hotspot_cnt-prev_hotspot_cnt, cur_hotspot_cnt))
        }
        
        cur_avg_geno_blocks_breaks[s] = 1
        prev_dist_cs = cur_avg_geno_dist_cumsum[s]
        prev_hotspot_cnt = cur_ns_rec_hotspot_cumsum$y[s]
        
      }
      else
      {
        #if (s>1){
        #  print(sprintf("%d), %f , %f",s, cur_dist_cs-prev_dist_cs, cur_hotspot_cnt-prev_hotspot_cnt))
        #}
      }  
    }
    
    diff(cur_avg_geno_dist_cumsum[cur_avg_geno_blocks_breaks>0.9])
    
    cur_avg_geno_blocks_breaks_tmp = cur_avg_geno_blocks_breaks;
    cur_avg_geno_blocks_breaks_tmp[1] = cur_avg_geno_blocks_breaks_tmp[1] + cur_max_black_num
    
    
    cur_avg_geno_blocks = cumsum(cur_avg_geno_blocks_breaks_tmp)
    
    cur_max_black_num = max(cur_avg_geno_blocks, na.rm = FALSE) 
    
    
    #cur_avg_geno_dist_cumsum
    
    block_avg_dist_cumsum = tapply(cur_avg_geno_dist_cumsum, as.factor(cur_avg_geno_blocks), mean )
    }
  
  cur_blocks_num = as.numeric(names(block_avg_dist_cumsum))
  cur_blocks_avg_cumsum = as.numeric(block_avg_dist_cumsum)
  cur_blocks_avg_dist_from_prev = c(Inf, diff(cur_blocks_avg_cumsum))
  
  cur_block_inds = match(cur_avg_geno_blocks,cur_blocks_num)
  
  avg_distance_FromPrevBlock = cur_blocks_avg_dist_from_prev[cur_block_inds]
  
  
  
  
  # adding to output
  geno_df$genetic_dist_log[cur_geno_inds] = log(cur_avg_geno_dist)
  geno_df$genetic_dist_cumsum[cur_geno_inds] = cur_avg_geno_dist_cumsum
  geno_df$genetic_rec_hotspot[cur_geno_inds] = cur_ns_rec_hotspot_num
  geno_df$genetic_block_breaks[cur_geno_inds] = cur_avg_geno_blocks_breaks
  geno_df$genetic_blocks[cur_geno_inds] = cur_avg_geno_blocks
  geno_df$genetic_avg_distance_FromPrevBlock[cur_geno_inds] = avg_distance_FromPrevBlock
  
}



# with header

write.table(geno_df,
            file = gzfile(output_filename),
            sep='\t' ,quote = FALSE, row.names = FALSE)














