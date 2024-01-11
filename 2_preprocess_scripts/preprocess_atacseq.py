
import argparse as ap

# merge peaks for atac-seq based data 
parser = ap.ArgumentParser(prog='This program processes in ATACseq data. Bins peaks into 600 bp bins')
parser.add_argument('-filepath', type=str, help='Input filepath of Excel file to process', nargs=1)
args = vars(parser.parse_args())