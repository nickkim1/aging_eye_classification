import argparse as ap
from Bio import SeqIO
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

parser = ap.ArgumentParser(prog='This program processes sequences generated from ATAC-seq data using bedtools getfasta')
parser.add_argument('-input_bed_file', type=str, help='Input bed file filepath used to generate fasta file', nargs=1)
parser.add_argument('-input_fasta_file', type=str, help='Input outputted fasta file filepath', nargs=1)
parser.add_argument('-output_csv_file', type=str, help='Input filepath of the csv file to write results to', nargs=1)
args = vars(parser.parse_args())

bed_file = args['input_bed_file'][0]
fasta_file = args['input_fasta_file'][0]
csv_file = args['output_csv_file'][0]

# read in bed file into a pandas dataframe. Read the counts into a numpy array 
bed = pd.read_table(bed_file, header=None, names=['chromosome', 'start', 'end', 'count'])
fasta = pd.read_table(fasta_file, header=None, names=['chromosome:interval', 'sequence'])

# exploration 
merged = pd.concat([fasta['sequence'], bed['count']], axis=1)
print('First three rows: ', merged.head(3))
median_count = merged['count'].median()
print('This is the median count: ', median_count)

# class metrics
print('Number of ones', (merged['count'] > median_count).astype('int').sum())
print('Number of zeros', merged.shape[0] - (merged['count'] > median_count).astype('int').sum())

# get boolean array discretizing counts in the array 
merged['count'] = (merged['count'] > median_count).astype('int')
print('Discretized df head', merged.head(3))

# add some code to concatenate the different files ** (ASSUMES WE HAVE MULTIPLE CELL TYPES OR SOMETHING SIMILAR) 

# write the resulting dataframes to CSV files 
train, test = train_test_split(merged, test_size=0.2)
train.to_csv(csv_file + '_train' + ".csv", index=False)
test.to_csv(csv_file + '_test' + ".csv", index=False)

