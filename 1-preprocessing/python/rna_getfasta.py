#  write the fasta file to a csv 
import argparse as ap
from Bio import SeqIO
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

parser = ap.ArgumentParser(prog='This program processes sequences generated from ATAC-seq data using bedtools getfasta')
parser.add_argument('-input_bed_file', type=str, help='Input bed file filepath used to generate fasta file', nargs=1)
parser.add_argument('-input_fasta_file', type=str, help='Input outputted fasta file filepath', nargs=1)
parser.add_argument('-output_csv_file', type=str, help='Input filepath of the csv file to write results to', nargs=1)
parser.add_argument('-output_csv_file_identifier', type=str, help='Input the identifier for your csv file to differentiate from other test/train files')
args = vars(parser.parse_args())

bed_file = args['input_bed_file'][0]
fasta_file = args['input_fasta_file'][0]
csv_file = args['output_csv_file'][0]
identifier = args['output_csv_file_identifier']

# read in bed file into a pandas dataframe. Read the counts into a numpy array 
bed = pd.read_table(bed_file, header=None, names=['chromosome', 'start', 'end', 'count'])
fasta = pd.read_table(fasta_file, header=None, names=['chromosome:interval', 'sequence'])

# exploration 
merged = pd.concat([fasta['sequence'], bed['count']], axis=1)
print('First three rows: ', merged.head(3))

# class metrics, already discretized 
print('Number of ones in whole dataset', (merged['count'] == 1).astype('int').sum())
print('umber of zeros in whole dataset', merged.shape[0] - (merged['count'] == 1).astype('int').sum())

# write the resulting dataframes to CSV files 
train, test = train_test_split(merged, test_size=0.2)
print('Number of ones in train dataset', (train['count'] == 1).astype('int').sum())
print('Number of zeros in train dataset', train.shape[0] - (train['count'] == 1).astype('int').sum())
print('Number of ones in test dataset', (test['count'] == 1).astype('int').sum())
print('Number of zeros in test dataset', test.shape[0] - (test['count'] == 1).astype('int').sum())
train.to_csv(csv_file + identifier + '_train' + ".csv", index=False, header=False)
test.to_csv(csv_file + identifier + '_test' + ".csv", index=False, header=False)

