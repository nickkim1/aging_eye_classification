import argparse as ap
import gffutils
import pyranges as pr
import pandas as pd
import numpy as np

# merge peaks for rnaseq data
parser = ap.ArgumentParser(prog='This program processes in RNAseq data into coverage format, binning into 600 bp bins')
parser.add_argument('-gtf_filepath', type=str, help='Input filepath of GTF file to process', nargs=1)
parser.add_argument('-counts_filepath', type=str, help='Input filepath of counts file to process', nargs=1)
parser.add_argument('-sequence_length', type=int, help="Input the length of the sequence you want for each gene")
args = vars(parser.parse_args())
gtf_filepath = args['gtf_filepath'][0]
counts_filepath = args['counts_filepath'][0]
sequence_length = args['sequence_length']

# first, load in the gtf file 
gr = pr.read_gtf(gtf_filepath)

# convert to pandas dataframe 
df = gr.df

# filter dataframe by "gene"
df = df[df['Feature'] == 'mRNA']
# print(df.head(5).iloc[:,4:]) 

# now retrieve the TSS coordinates. 
    # * start/end coords are given for DNA strand 
    # * 1 - if gene is on + (coding) strand, then it is oriented 5' -> 3' and so start position on xsome aligns w/ TSS of gene (which is transcribed 5'->3')
    # * 2 - if gene is on - (template) strand, then it is oriented 3' -> 5' so end position on chromosome is on the 5' end of gene and is actually where transcription starts
    # thus, + has TSS = start. - has TSS = end.  

tss_sites = pd.Series([df.iloc[i, 3] if df.iloc[i, 6] == "-" else df.iloc[i, 4] for i in range(df.shape[0])])

# setup new start and end coordinates, set this to whatever you want 
tss_relative_starts = pd.Series(tss_sites - (sequence_length // 2))
tss_relative_ends = pd.Series(tss_sites + (sequence_length // 2))
tss_relative_starts = tss_relative_starts.astype(int)
tss_relative_ends = tss_relative_ends.astype(int)

# now create a new dataframe with the chromosome, gene id, start and stop 
tss_df = pd.DataFrame({
    "Chromosome": df["Chromosome"].values,
    "gene_name": df["gene_id"].values, 
    "Start": tss_relative_starts, 
    "Stop": tss_relative_ends
})

# protocol: 
    # we generically generated tss starts/ends from every gene in the D. melanogaster genome
    # now, we need to check for matching gene ID across the generic DF + each DF from each timepoint in excel file to assign the right counts

# now, we need to read in the counts matrix and just find matching gene names
# the principle is that because we read in the GTF, we have a reference dataset that we can essentially use to map TSS start/end bounds to counts

counts_df = pd.read_csv(filepath_or_buffer=counts_filepath, sep='\t')

# we want to average across all the replicates within the dataframe. select everything excluding the first column 
avged_rep_counts = pd.Series(counts_df.iloc[:, 1:].mean(axis=1))
counts_df = counts_df.assign(avged_replicate_counts=avged_rep_counts)

# discretize w.r.t median count within the avgg replicates column first 
median_count = counts_df.iloc[:, counts_df.shape[1]-1].median(axis=0)
binarized_counts = pd.Series(np.where(avged_rep_counts > median_count, 1, 0))
binarized_counts = binarized_counts.astype(int)
counts_df = counts_df.assign(binarized_median_counts=binarized_counts)
counts_df = counts_df.filter(items=['gene_name', 'binarized_median_counts'])

# now merge with the gtf dataframe 
merged_tss_counts_df = pd.merge(tss_df, counts_df, on='gene_name', how='inner')
merged_tss_counts_df = merged_tss_counts_df.filter(items=["Chromosome", "Start", "Stop", "binarized_median_counts"])

# adjust the following line, but the idea is to write the above dataframe to a bed file.
merged_tss_counts_df.to_csv(f'../0-data/1-experimental/bed/agingeye/cts_bed.bed', sep='\t', index=False, header=False)


# old printlines: 
# print(">>>> Shape of the dataframe:", df.shape, "\n")
# print(">>>> Column ids of the dataframe: \n", df.columns, "\n")
# print(">>>> First row of the dataframe: \n", df.iloc[0,:])
# print(">>>> Unique feature identifiers in PD dataframe:", pd.unique(df['Feature']), "\n")
# print(">>>> Number of mRNA transcripts within the GTF file:", len(df[df['Feature'] == 'mRNA']))
# do some feature exploration 
# print("These are the first five entries of the GTF file, in dataframe format")
# print(df.head(5))
# print("These are now the first five entries of the GTF file, in dataframe format, after filtering for only transcripts (only show last couple of columns to show gnes)")
# print("This is a sample of the appropriate TSS sites in the dataframe")
# print(tss_sites.head(5))
# print("This is the original df post-TSS concatenation")
# print(tss_df.head(2)) 
# print("These are its dimensions:", tss_df.shape)
# print("This is the shape of the counts df: ", counts_df.shape)
# print(counts_df.head(6), counts_df.shape)
# print("This is the median count", median_count)
# print("Verification of median count", avged_rep_counts.median())
# print(merged_tss_counts_df.head(2))