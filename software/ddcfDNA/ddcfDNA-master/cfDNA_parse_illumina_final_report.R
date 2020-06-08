

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
# setting input parameters
######################################################################

# setting the value of the arguments
in_final_report_filename <- args[1]
patient_id <- args[2]
in_manifest_filename <- args[3]
out_genotype_filename <- args[4]

#in_final_report_filename <- '/workdir/apc88/ddcfDNA/Data/Vlaminck_RQ013504_FINAL_101918_FinalReport.txt'
#patient_id <- 'R2'
#in_manifest_filename <- '/workdir/apc88/ddcfDNA/workflow/Input/illumina/HumanOmni25.tsv'
#out_genotype_filename <- 'Output/genotypes/2.genotype.tsv.gz'

print(sprintf("in_final_report_filename:%s",in_final_report_filename))
print(sprintf("patient_id:%s",patient_id))
print(sprintf("in_manifest_filename:%s",in_manifest_filename))
print(sprintf("out_genotype_filename:%s",out_genotype_filename))

######################################################################
# parsing params
######################################################################

#patient_id_str = as.character(as.integer(patient_id))
#patient_id = as.integer(patient_id)
patient_id_str=patient_id
patient_id=patient_id
######################################################################
# loading librares
######################################################################

# loading libraries


if (!require("reshape2")) {
  library(reshape2)
}


######################################################################
# util functions
######################################################################

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

flip_nt <- function(nt){


  if (nt == 'A' | nt == 'a') {
    ret_nt = 'T'
  } else if (nt == 'G' | nt == 'g') {
    ret_nt = 'C'
  } else if (nt == 'C' | nt == 'c') {
    ret_nt = 'G'
  } else if (nt == 'T' | nt == 't') {
    ret_nt = 'A'
  } else {
    print(sprintf(" Can NOT flip: %s", nt))
    ret_nt = nt
  }

  ret_nt
}

filter_nt <- function(nt_vec1,nt_vec2){

  nt_list = c("A","G","C","T")

  ret_filter = logical(length=length(nt_vec1))

  for (i in 1:length(nt_vec1)) {
    cur_f = (nt_vec1[i] %in% nt_list) & (nt_vec2[i] %in% nt_list)
    if (!cur_f) {
      print(sprintf("NOT ok nts, line %d: %s,%s\n", i, nt_vec1[i],nt_vec2[i]))
    }
    ret_filter[i] = cur_f
    if (i %% 10000 == 0) {
      print(sprintf("filter nt line %d\n", i))
    }
  }
  ret_filter
}


######################################################################
# extracting the patient data from the full report table
######################################################################


# saving to a tmp file without the header row

patient_full_Report_tabel_namefile = paste(in_final_report_filename,".",patient_id_str,sep="")
in_fh = file(in_final_report_filename,"rt")
out_fh = file(patient_full_Report_tabel_namefile, "w")

cur_l= 0
data_header = "[Data]"
max_lines_for_header = 20
is_data_flag = FALSE
for (l in 1:max_lines_for_header) {
  input <- trim(readLines(in_fh, n=1))
  cur_l=cur_l+1
  print(sprintf("Header line:%s", input))
  if (input == data_header) {
    is_data_flag=TRUE
    break
  }
}

if (!is_data_flag) {
  write("ERROR full report - data header line does not appear in first 20 lines", stderr())
  stop("ERROR full report - data header line does not appear in first 20 lines")
} else {

  w_lines = 1
  cur_l = cur_l + 1
  # wrting the header
  writeLines(readLines(in_fh, n=1), con=out_fh)

  while (length(input<- readLines(in_fh, n=65536))>  0){

    cur_l = cur_l + length(input)

    tmp = strsplit(input,"\t")
    tmp_ind = sapply(tmp, function(x) x[2]==patient_id_str)

    writeLines(input[tmp_ind], con=out_fh)

    w_lines = w_lines +sum(tmp_ind)

    print(sprintf("full report read %d lines, wrote %d lines",cur_l, w_lines))
  }
}

close(in_fh)
close(out_fh)

#################################################################
# loading tables, merging and parsing
#################################################################

final_report = read.delim(patient_full_Report_tabel_namefile)
manifest = read.delim(in_manifest_filename)

colnames(final_report)
colnames(manifest)

if ("Chr" %in% colnames(final_report) &&
    "Position" %in% colnames(final_report)){
  # modiy final report chr name to match the manifest and fasta file
  fr_chr_vec = as.vector(final_report$Chr)
  fr_chr_vec[fr_chr_vec=="MT"] = "M"
  fr_chr_vec = paste("chr",fr_chr_vec, sep="")


  final_report[, 'Chr'] = factor(fr_chr_vec)


  merged_table = merge(final_report,manifest, by=c("Chr","Position") )

} else if ("SNP.Name" %in% colnames(final_report))
{
  if ("Chr" %in% colnames(final_report))
  {
    final_report$Chr = NULL
  }
  if ("Position" %in% colnames(final_report))
  {
    final_report$Position = NULL
  }

  final_report$merge_snp_id = final_report$SNP.Name
  manifest$merge_snp_id = manifest$Name

  merged_table = merge(final_report,manifest, by=c("merge_snp_id") )

} else {
  stop(sprintf("Can not join manifest and final report tables #### %s \n##### %s",
               paste(colnames(manifest), collapse = ' '),
               paste(colnames(final_report), collapse = ' ') )  )

}


allele1_top_colname = "Allele1...Top"
allele2_top_colname = "Allele2...Top"
gc_score_colname = "GC.Score"
flip_top_colname = "is_flip_top"

print(sprintf("Number of row final report %d\n", nrow(final_report)) )
print(sprintf("Number of row manifest %d\n", nrow(manifest)) )
print(sprintf("Number of row after joining %d\n", nrow(merged_table)) )

# flip the top alleles

nt_allele1 = merged_table[,allele1_top_colname]
nt_allele2 = merged_table[,allele2_top_colname]
filp_top = merged_table[,flip_top_colname]==1

# filt indices
ok_nts = c("A","G","C","T")
row_filt = (nt_allele1 %in% ok_nts) & (nt_allele2 %in% ok_nts)

# flipping the nt in minus strand
nt_allele1_flipped = nt_allele1
nt_allele2_flipped = nt_allele2

cur_inds = filp_top & nt_allele1 == "A"
nt_allele1_flipped[cur_inds] = "T"
cur_inds = filp_top & nt_allele1 == "G"
nt_allele1_flipped[cur_inds] = "C"
cur_inds = filp_top & nt_allele1 == "C"
nt_allele1_flipped[cur_inds] = "G"
cur_inds = filp_top & nt_allele1 == "T"
nt_allele1_flipped[cur_inds] = "A"

cur_inds = filp_top & nt_allele2 == "A"
nt_allele2_flipped[cur_inds] = "T"
cur_inds = filp_top & nt_allele2 == "G"
nt_allele2_flipped[cur_inds] = "C"
cur_inds = filp_top & nt_allele2 == "C"
nt_allele2_flipped[cur_inds] = "G"
cur_inds = filp_top & nt_allele2 == "T"
nt_allele2_flipped[cur_inds] = "A"


print(sprintf("Total SNPs %d, ok SNPs %d, not ok %d, need flipping %d",length(row_filt),sum(row_filt),length(row_filt)-sum(row_filt), sum(filp_top) ))


out_df = data.frame(
  Chr = factor(merged_table$Chr[row_filt]),
  Position = factor(merged_table$Position[row_filt]),
  Allele1 = factor(nt_allele1_flipped[row_filt]),
  Allele2 = factor(nt_allele2_flipped[row_filt]),
  GCScore = factor(merged_table[row_filt,gc_score_colname]) )

write.table(out_df,
          file = gzfile(out_genotype_filename),
          sep='\t' ,quote = FALSE, row.names = FALSE)

print(sprintf("wrote output to:%s", out_genotype_filename))
