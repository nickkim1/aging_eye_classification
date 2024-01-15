# --- script for getting list of TSS regions from a provided gtf file --- # 

# set directorie
setwd("/oscar/data/rsingh47/nhkim/IISAGE/Nick-Basset")

# install packages, restart R session (if necessary) for everything to get imported correctly w/o errors 
BiocManager::install("rtracklayer") # i don't have this currently installed in the container i'm working out of
library(rtracklayer)
library(dplyr)
library(ggplot2)
install.packages('seqinr')
library(seqinr)
library("Biostrings")

# load in gtf as a granges object, convert to df 
gtf <- rtracklayer::import('./preprocessing/data/refs/Drosophila_melanogaster.BDGP6.46.110.gtf')
gtf_df=as.data.frame(gtf) 

# also load in the counts matrix 
counts_df <- read.table('./preprocessing/data/aging_eye/GSE202053_RNAseq_counts.txt', header=TRUE)

# filter for ONLY genes in the gtf -> get tss sites of those
tss_sites <- dplyr::filter(gtf_df, gtf_df$type=='gene')
tss_sites$tss <- ifelse(tss_sites$strand == "+", tss_sites$start, tss_sites$end) # IF positive -> tss = start. IF negative -> tss = end. Adds new tss column to the df.

# create a temporary bed file based on this tss information 
tss_bed_df <- data.frame(chrom=tss_sites$seqnames, gene_name=tss_sites$gene_id, chromStart=tss_sites$start-300, chromEnd=tss_sites$start+300)
merged_tss_df <- merge(counts_df, tss_bed_df, by="gene_name") 

# ---------------------- TESTING IGNORE ---------------------- # 
bed_df_gn <- tss_bed_df$gene_name
counts_gn <- counts_df$gene_name
common_gn <- intersect(bed_df_gn, counts_gn)
# Find elements in df1$ColumnA not present in df2$ColumnA
not_in_df2 <- counts_df$gene_name[!(counts_df$gene_name %in% tss_bed_df$gene_name)]
print(not_in_df2)
# ------------------------------------------------------------ # 

non_zero_indices <- c(as.numeric(which(merged_tss_df[,"chromStart"]<0)))
for (index in non_zero_indices) {
  print(index)
  merged_tss_df[index,"chromStart"] <- 0
  merged_tss_df[index, "chromEnd"] <- 600
}

# double check -- EXCLUSIVELY for the aging-eye data 
# print(merged_tss_df[1305,])

# select for only the required fields in the bed file 
merged_tss_df <- dplyr::select(merged_tss_df, chrom, chromStart, chromEnd)

# write to preprocessing directory (NOTE: 53 genes are excluded from the counts matrix...apparently not in the gtf?? no clue why)
write.table(x=merged_tss_df,file="./preprocessing/data/final_tss_bed.bed", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)

# now read in the fasta file with the sequence information using BioStrings -> dataframe 
preprocessed_tss_fasta <- readDNAStringSet("./preprocessing/data/final_bed.fa")
seq_name = names(preprocessed_tss_fasta)
sequence = paste(preprocessed_tss_fasta)
preprocessed_fasta_df <- data.frame(seq_name, sequence)

# now we want to merge gene names -> preprocessed_fasta_df, with corresponding counts for each condition
in_counts <- as.vector(which(counts_df$gene_name %in% tss_bed_df$gene_name))
# we know there are 53 entries excluded (at least for this data), so merge based on indices common to counts df
# this works bc preprocessed_fasta_df is of the same length as the original, merged, df w/ 1:1 correspondence
subsetted_counts <- counts_df[in_counts,]
preprocessed_fasta_df <- cbind(gene_name=subsetted_counts$gene_name, sequence=preprocessed_fasta_df[,2:ncol(preprocessed_fasta_df)], Set2R1=subsetted_counts$Set2R1, Set2R2=subsetted_counts$Set2R2, Set2R3=subsetted_counts$Set2R3)

# write to a new file for the new counts matrix 
write.csv(x=preprocessed_fasta_df, file="./preprocessing/data/final_counts.csv", quote=FALSE, row.names = FALSE)

# ----------------------- OLD CODE IGNORE ----------------------------------- # 
# TODO: 
# -- idea is now to run getfasta on this -> extract out corresponding subsequence, put that in a dataframe
# -- for each entry in the dataframe, then just put in corres. count (labels) -> write to table -> input to model
 
# subset_df <- temp_bed_file %>% slice(2:4)

# first, write the table to a bed file
# testing <- data.frame(c("3r", 100, 200))
# write.table(x=temp_bed_file,file="./preprocessing/data/final_get_fasta_bed.bed", quote=FALSE, sep="\t", col.names = FALSE, row.names = FALSE)
# write.table(x=testing,file="./preprocessing/data/test3.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




