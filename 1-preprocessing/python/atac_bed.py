import pybedtools as pb
import argparse as ap
import os 
import numpy as np

# Credits for conceptualization of the algorithm go to the authors of the original Basset paper. 
# This is mostly just a modification of what they had already implemented. 

# merge peaks for atac-seq based data 
parser = ap.ArgumentParser(prog='This program processes in ATACseq data in BED format. Bins peaks into 600 bp bins')
parser.add_argument('-input_bed_directory', type=str, help='Input directory holding BED files to process, in comma delimited format', nargs=1)
parser.add_argument('-output_bed_directory', type=str, help='Input directory to hold processed BED files, in comma delimited format')
parser.add_argument('-max_overlap', type=int, help='Maximum overlap between consecutive intervals in the BED file')
parser.add_argument('-sequence_length', type=int, help="How long you want input sequences to the model to be")
args = vars(parser.parse_args())

# set global input arguments to be used in the script
directory = args['input_bed_directory'][0]
output_directory = args['output_bed_directory']
max_overlap = args['max_overlap']
sequence_length = args['sequence_length']

# ------------------------------------------------------------------------------ # 
# Methods for generating merged results 
#------------------------------------------------------------------------------ # 

# generate active input metrics 
def generate_metrics(file_dict:dict, chroms_starts_ends:dict, file, unique_intervals:dict):
    for row in file: 
        file_dict[row.chrom] = file_dict.get(row.chrom, {})
        file_dict[row.chrom][(int(row[1]), int(row[2]))] = [float(row[4])]
        unique_intervals[(int(row[1]), int(row[2]))] = 1 + unique_intervals.get((int(row[1]), int(row[2])), 0)
    for chromosome, all_rows in file_dict.items(): 
        # key: chromosome -> (start, stop)
        chroms_starts_ends[chromosome] = (min(list(all_rows), key = lambda x: x[0])[0],  max(list(all_rows), key = lambda x : x[1])[1])
    return file_dict

# ------------------------------------------------------------------------------- # 
# Class definition for Peak data 
# ------------------------------------------------------------------------------- # 
class Peak():
    def __init__(self, start, end, count, weight): 
       self.start = start
       self.end = end
       self.mid = self.calculate_mid(start, end)
       self.count = int(count)
       self.weight = weight
    
    def extend_coords(self, mid, sequence_length, chromosome, chroms_starts_ends): 
        self.start = int(mid - sequence_length / 2)
        self.end = int(mid + sequence_length / 2)
        # check upper and lower bounds 
        if self.start < chroms_starts_ends[chromosome][0]: 
            self.start, self.end = chroms_starts_ends[chromosome][0], chroms_starts_ends[chromosome][0] + 600
        # this condition essentially means i'll probably just merge the last peak later on 
        if self.end > chroms_starts_ends[chromosome][1]: 
            self.start, self.end = chroms_starts_ends[chromosome][1] - 600, chroms_starts_ends[chromosome][1] 

    def calculate_mid(self, s, e):
        return s + ((e - s) / 2)
    
    def update_counts(self, a, b):
        self.count = int(a + b)

    # function that merges any given new peak with the specified previous one 
    def merge(self, peak_to_merge_with, sequence_length, chromosome, chroms_starts_ends): 
        
        cur_interval, interval_to_merge = (self.start, self.end), (peak_to_merge_with.start, peak_to_merge_with.end)
        cur_mid, prev_mid = self.calculate_mid(cur_interval[0], cur_interval[1]), self.calculate_mid(interval_to_merge[0], interval_to_merge[1])
        cur_weight, prev_weight = self.weight, peak_to_merge_with.weight 
        new_merge_midpt = (cur_weight * cur_mid + prev_weight * prev_mid) / (cur_weight + prev_weight)

        # update counts, coords
        self.update_counts(self.count, peak_to_merge_with.count)
        self.extend_coords(new_merge_midpt, sequence_length, chromosome, chroms_starts_ends)


# merge the peaks, in each chromosome, so that they are 600 bp long. followed procedure from the paper 
def extend_peaks(file_dict, row_dict, chroms_starts_ends, sequence_length):
    
    # 0. for each chromosome => 
    # 1. iterate thru peaks, extending each one to +/- y bp from midpoint. 

    for chromosome, all_rows in file_dict.items():

        # sort all intervals by start coordinate in case they weren't already 
        sorted_intervals = sorted(all_rows.keys(), key=lambda x : x[0])
        # create a peaks list for all the intervals in the sorted intervals list 
        sorted_peaks_list = [Peak(interval[0], interval[1], all_rows[(interval[0], interval[1])][0], int(all_rows[(interval[0], interval[1])][1])) for interval in sorted_intervals]
        # iterate over all the intervals until you encounter a peak that overlaps > 200 bp with the peak you're currently at. the objective is to merge that into our current peak
        for i, cur_peak in enumerate(sorted_peaks_list):
            # for every peak just extend it. update entry in sorted_peaks_list 
            cur_peak.extend_coords(cur_peak.mid, sequence_length, chromosome, chroms_starts_ends) 

        # TESTING ~~~ (check if all sequences are 600 bp long) 
        var = [(peak.start, peak.end) for peak in sorted_peaks_list if peak.end - peak.start != 600]
        if len(var) > 0: 
            print(var)
        
        # update the row_dict w/ the intervals that are newly merged. now can iterate thru the dict to create data!  
        row_dict[chromosome] = sorted_peaks_list

    return row_dict

def greedily_merge_peaks(max_overlap, row_dict, chroms_starts_ends, sequence_length):

    # => now greedily merge the extended sequences to not be within the specified overlap distance of each other
    # key -> chromosome: value -> sorted list of all peaks that belong to that chromosome 
    for chromosome in row_dict.keys():
        all_peaks = row_dict[chromosome]
        i = 0
        while i < len(all_peaks)-1:  
            # check if it overlaps with next peak, if so, merge the two peaks and check again 
            if all_peaks[i].end - all_peaks[i+1].start > max_overlap: 
                # merge next peak, so that when u iterate to next peak ur evaluating the new peak's overlap with the i+2 one
                all_peaks[i+1].merge(peak_to_merge_with=all_peaks[i], sequence_length=sequence_length, chromosome=chromosome, chroms_starts_ends=chroms_starts_ends)
                all_peaks = all_peaks[:i] + all_peaks[i+1:]
            i += 1
    return row_dict 

# notes how active each individual peak, based on its presence across all cell types. for me, age points? 
def input_active_peaks(file_dict:dict, unique_intervals:dict):

    # if we're testing by pure equality alone, then just use the existing dict
    for rows_per_chromosome in file_dict.values(): 
        # should just be 1 interval per assuming non-overlapping windows 
        for interval in rows_per_chromosome.keys():
            rows_per_chromosome[interval].append(unique_intervals.get(interval, 1))

    return file_dict

def write_dict_to_bed_format(row_dict:dict, output_directory, filepath_name):
    concat_path = output_directory + filepath_name
    with open(concat_path, 'w') as bed_file: 
        for chrom, rows in row_dict.items(): 
            for peak in rows: 
                chromStart = str(peak.start)
                chromEnd = str(peak.end)
                readCount = str(peak.count)
                bed_line = "\t".join([chrom, chromStart, chromEnd, readCount])
                bed_file.write(bed_line + '\n')   

def main(directory, output_directory, sequence_length):
    if os.path.exists(directory) and os.path.isdir(directory):
        # iterate thru all files 
        unique_intervals = {}
        all_files = []
        for file in os.listdir(directory):
            filepath = os.path.join(directory, file) 
            if os.path.isfile(filepath):
                # generate dictionary for the file features
                chroms_starts_ends = {}
                file_dict = generate_metrics({}, chroms_starts_ends, pb.BedTool(filepath), unique_intervals)
                all_files.append((file, file_dict))
        # iterate thru files again to update the unique interval counts, greedily merge peaks, write output in BED format to new file 
        for (file, file_dict) in all_files: 
            file_dict = input_active_peaks(file_dict, unique_intervals)
            row_dict = extend_peaks(file_dict, {}, chroms_starts_ends, sequence_length)
            row_dict = greedily_merge_peaks(max_overlap, row_dict, chroms_starts_ends, sequence_length)
            write_dict_to_bed_format(row_dict, output_directory, file)
            print("File ", file, " processed and outputted to new BED file!")
    else: 
        print('Inputted directory filepath does not exist')

if __name__ == '__main__':
    main(directory, output_directory, sequence_length)


# ---------------------------- IGNORE ---------------------------- # 
# else: 
#     # max overlap isn't exceeded, so just append the peaks to this list. you'll have a separate method that iterates through this list and resets it once you're done with it
#     largest_end = max(largest_end, interval[1])
#     peaks_to_merge.append(())
    
# # get all intervals btwn peak before extension and after extension
# largest_end = max(largest_end, sorted_intervals[i][1])
# # print(cur_peak.end, largest_end)
# peaks_to_merge.append(sorted_intervals[i])
# i += 1 
    
# TESTING ~~~
# if sorted_intervals[i][0] == 20801140:
#     print("Original coords: ", sorted_intervals[i][0], sorted_intervals[i][1])
#     # print("Sorted peaks list", [print(peak) for peak in sorted_peaks_list if peak.start == 20801140])
#     for k in range(len(sorted_peaks_list)):
#         if sorted_peaks_list[k].start == 20801140: 
#             print(k)
#             # print(sorted_peaks_list[k-1].start, sorted_peaks_list[k-1].end)
#     for m in range(len(sorted_intervals)):
#         if sorted_intervals[m][0] == 20801140:
#             print(m)
#             # print(sorted_intervals[k-1][0], sorted_intervals[k-1][1])
#     print("Extended coords: ", cur_peak.start, cur_peak.mid, cur_peak.end)
#     print("Previous peak", sorted_peaks_list[i-1].start, sorted_peaks_list[i-1].end)

# print("Coord i: ", peaks_to_merge[i][1].start, peaks_to_merge[i][1].end)
# print("Coord i+1: ", peaks_to_merge[i+1][1].start, peaks_to_merge[i+1][1].end)
    
# print(peaks_to_merge[i+1][0])

# while len(peaks_to_merge) > 1 and cur_max_overlap >= max_overlap: 

# take into account your current peak, merge based on midpoints
# cur_max_overlap = max_overlap

# --- OLD 
# # the condition before guaranteed overlap <= 200 bp using UNMERGED peaks. 
# # but after you merge peaks, you will likely have additional overlap because you're extending each peak +/- 300 bp

# # then, if the max_overlap you find is >= specified max_overlap then merge. true for first peak always once inputted, but maybe not true for peaks post-merging 
# if cur_max_overlap >= max_overlap:
#     # first, try to find max_overlap in peaks to merge. in the beginning, the max_overlap will be guaranteed <= cur_max_overlap. but as you keep merging that might not be the case and you might not need to merge anymore
#     merge_pter = 0
#     for i in range(1, len(peaks_to_merge)-1): 
#         # if negative, that means NO overlap. if positive, then YES overlap. hence greater than sign 
#         temp_overlap = peaks_to_merge[i].end - peaks_to_merge[i+1].start
#         if temp_overlap > cur_max_overlap: 
#             cur_max_overlap = temp_overlap
#             merge_pter = i
 # new midpt is: weighted avg of midpts of the peaks. get new params for the new peak
# cur_interval, interval_to_merge = (peaks_to_merge[merge_pter].start, peaks_to_merge[merge_pter].end), (peaks_to_merge[merge_pter+1].start, peaks_to_merge[merge_pter+1].end)
# cur_mid, next_mid = self.calculate_mid(cur_interval[0], cur_interval[1]), self.calculate_mid(interval_to_merge[0], interval_to_merge[1])
# a, b = peaks_to_merge[merge_pter].weight, peaks_to_merge[merge_pter+1].weight
# mid = (a * cur_mid + b * next_mid) / (a + b)

# # update the new merged peak. 
# peaks_to_merge[merge_pter].update_counts(peaks_to_merge[merge_pter].count, peaks_to_merge[merge_pter+1].count)
# peaks_to_merge[merge_pter].extend_coords(mid, sequence_length)

# # ~ iterate by popping off list for the peak you just incorporated 
# peaks_to_merge = peaks_to_merge[:merge_pter+1] + peaks_to_merge[merge_pter+2:]

# # ~ keep iterating. Should modify sorted_intervals in place 

# # new midpt is: weighted avg of midpts of the peaks. get new params for the new peak
# cur_interval, interval_to_merge = (peaks_to_merge[merge_pter].start, peaks_to_merge[merge_pter].end), (peaks_to_merge[merge_pter+1].start, peaks_to_merge[merge_pter+1].end)
# cur_mid, next_mid = self.calculate_mid(cur_interval[0], cur_interval[1]), self.calculate_mid(interval_to_merge[0], interval_to_merge[1])
# a, b = peaks_to_merge[merge_pter].weight, peaks_to_merge[merge_pter+1].weight
# mid = (a * cur_mid + b * next_mid) / (a + b)

# # update the new merged peak. 
# peaks_to_merge[merge_pter].update_counts(peaks_to_merge[merge_pter].count, peaks_to_merge[merge_pter+1].count)
# peaks_to_merge[merge_pter].extend_coords(mid, sequence_length)

# # ~ iterate by popping off list for the peak you just incorporated 
# peaks_to_merge = peaks_to_merge[:merge_pter+1] + peaks_to_merge[merge_pter+2:]

# # ~ keep iterating. Should modify sorted_intervals in place 