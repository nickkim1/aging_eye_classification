import csv 
import argparse as ap 
import pandas as pd

# parse in command line inputs 
parser = ap.ArgumentParser(prog='This program processes in DESeq RNAseq data in Excel format')
parser.add_argument('-filepath', type=str, help='Input filepath of Excel file to process', nargs=1)
parser.add_argument('-sheetnames', type=str, help='Input names of sheet names in the Excel file in comma-delimited format')
parser.add_argument('-skiprows', type=str, help='Input number of rows to skip in comma-delimited format')
parser.add_argument('-header', type=str, help='Input index of the header to use in the CSV file')
args = vars(parser.parse_args())

# parse in argument values 
filepath = args['filepath'][0]
sheet_name = args['sheetnames'].split(',')
# NOTE: skip rows, header assumes a standardized # of skipped rows, header idx across all workbooks 
skip_rows = args['skiprows'].split(',')
skip_rows = [int(s) for s in skip_rows]
print('=====================================')
print('These are the fields inputted: ', filepath, sheet_name, skip_rows)
print('=====================================')

# read in the Excel file into a pandas dataframe 
excel_file = pd.read_excel(io=filepath, sheet_name=sheet_name, skiprows=skip_rows)
for sheet_name, file in excel_file.items(): 
    print(f'2 row head of file: {sheet_name}')
    print(file.head(2))
    print('=====================================')

# the above can then be used to label the data 
    
