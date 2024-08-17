import argparse as ap
import gffutils
import pyranges as pr
import pandas as pd

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
print("=====================================================================================================================")
print(">>>> Shape of the dataframe:", df.shape, "\n")
print(">>>> Column ids of the dataframe: \n", df.columns, "\n")
print(">>>> First row of the dataframe: \n", df.iloc[0,:])
print(">>>> Unique feature identifiers in PD dataframe:", pd.unique(df['Feature']), "\n")
print(">>>> Number of mRNA transcripts within the GTF file:", len(df[df['Feature'] == 'mRNA']))

# do some feature exploration 
print("=====================================================================================================================")
print("These are the first five entries of the GTF file, in dataframe format")
print(df.head(5))

# filter dataframe by "gene"
print("=====================================================================================================================")
print("These are now the first five entries of the GTF file, in dataframe format, after filtering for only transcripts (only show last couple of columns to show gnes)")
df = df[df['Feature'] == 'mRNA']
print(df.head(5).iloc[:,4:]) 

# now retrieve the TSS coordinates. 
    # * start/end coords are given for DNA strand 
    # * 1 - if gene is on + (coding) strand, then it is oriented 5' -> 3' and so start position on xsome aligns w/ TSS of gene (which is transcribed 5'->3')
    # * 2 - if gene is on - (template) strand, then it is oriented 3' -> 5' so end position on chromosome is on the 5' end of gene and is actually where transcription starts
    # thus, + has TSS = start. - has TSS = end. 

tss_sites = pd.Series([df.iloc[i, 3] if df.iloc[i, 6] == "-" else df.iloc[i, 4] for i in range(df.shape[0])])
print("=====================================================================================================================")
print("This is a sample of the appropriate TSS sites in the dataframe")
print(tss_sites.head(5))

# setup new start and end coordinates, set this to whatever you want 
tss_relative_starts = tss_sites - (sequence_length // 2)
tss_relative_ends = tss_sites + (sequence_length // 2)

# now create a new dataframe with the chromosome, gene id, start and stop 
tss_df = pd.DataFrame({
    "Chromosome": df["Chromosome"].values,
    "gene_id": df["gene_id"].values, 
    "Start": tss_relative_starts, 
    "Stop": tss_relative_ends
})

print("=====================================================================================================================")
print("This is the original df post-TSS concatenation")
print(tss_df.head(5)) 
print("These are its dimensions:", tss_df.shape)

# protocol: 
    # we generically generated tss starts/ends from every gene in the D. melanogaster genome
    # now, we need to check for matching gene ID across the generic DF + each DF from each timepoint in excel file to assign the right counts

# TODO: hard encoded the sheet names, but should figure out a way to just pass this as a CLI argument 
sheet_names = ['D20vsD10', 'D30vsD10', "D40vsD10", "D50vsD10", "D60vsD10"]

# makes a dictionary of each sheet in the excel file, converted to a dataframe 
dfs = {sheet_name: pd.read_excel(counts_filepath, sheet_name, skiprows=2, index_col=None, na_values=["NA"])
          for sheet_name in sheet_names}

# iterate over the collection of dataframes (dfs)
for i, sheet_name in enumerate(dfs.keys()): 

    # merge counts by gene_ids 
    merged_result = pd.merge(dfs[sheet_name], tss_df, on="gene_id")

    # tenatively binarize counts, using LFC 1 as a threshold divider 
    merged_result["Count"] = (merged_result['log2FoldChange'] > 1.5).astype(int)

    # get rid of extraneous columns, just keep necessary ones for the bed file 
    adjusted_df = pd.DataFrame({
        "chrom": merged_result["Chromosome"].values, 
        "chromStart": merged_result["Start"].values,
        # "gene": merged_result["gene_id"].values, 
        "chromEnd": merged_result["Stop"].values, 
        "count": merged_result["Count"].values
    })

    print(f"These are the dimensions for the {i}'th sheet in the file", adjusted_df.shape)

    # write that dataframe to a txt file in bed format
    adjusted_df.to_csv(f'../0-data/1-experimental/bed/agingeye/{sheet_name}.bed', sep='\t', index=False, header=False)

# next, run getfasta to get the appropriate sequence for the intervals identified above 
