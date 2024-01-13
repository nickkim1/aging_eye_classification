import pybedtools as pb
import argparse as ap
import os 
import numpy as np

# merge peaks for atac-seq based data 
parser = ap.ArgumentParser(prog='This program processes in ATACseq data in BED format. Bins peaks into 600 bp bins')
parser.add_argument('-directory', type=str, help='Input directory holding BED files to process, in comma delimited format', nargs=1)
parser.add_argument('-maxoverlap', type=int, help='Maximum overlap between consecutive intervals in the BED file')
parser.add_argument('-sequencelength', type=int, help="How long you want input sequences to the model to be")
args = vars(parser.parse_args())

# set global input arguments to be used in the script
directory = args['directory'][0]
max_overlap = args['maxoverlap']
sequence_length = args['sequencelength']

# ------------------------------------------------------------------------------- # 
# Methods for generating merged results 
# ------------------------------------------------------------------------------- # 

# generate active input metrics 
def generate_metrics(file_dict:dict, file, unique_intervals:dict):
    for row in file: 
        file_dict[row.chrom] = file_dict.get(row.chrom, {})
        file_dict[row.chrom][(int(row[1]), int(row[2]))] = int(float(row[4]))
        unique_intervals[(int(row[1]), int(row[2]))] = 1 + unique_intervals.get((int(row[1]), int(row[2])), 0)
    return file_dict


# ------------------------------------------------------------------------------- # 
# Class definition for Peak data 
# ------------------------------------------------------------------------------- # 
class Peak():
    def __init__(self, start, end, count): 
       self.start = start
       self.end = end
       self.mid = self.calculate_mid(start, end)
       self.count = count
    
    def extend_coords(self, mid, sequence_length): 
        self.start = mid - sequence_length / 2
        self.end = mid + sequence_length / 2

    def calculate_mid(self, s, e):
        return s + ((e - s) / 2)
    
    def update_counts(self, a, b):
        self.count = a + b

    # function that merges any given new peak
    def merge(self, peaks_to_merge:list, max_overlap, sorted_peaks_list, sequence_length): 

        # take into account your current peak, merge based on midpoints
        cur_max_overlap = max_overlap

        # the condition before guaranteed overlap <= 200 bp using UNMERGED peaks. 
        # but after you merge peaks, you will likely have additional overlap because you're extending each peak +/- 300 bp

        while len(peaks_to_merge) > 1 and cur_max_overlap >= max_overlap: 

            # first, try to find max_overlap in peaks to merge. in the beginning, the max_overlap will be guaranteed <= cur_max_overlap. but as you keep merging that might not be the case and you might not need to merge anymore
            merge_pter = 0
            for i in range(len(sorted_peaks_list)-1): 
                # if negative, that means NO overlap. if positive, then YES overlap. hence greater than sign 
                temp_overlap = sorted_peaks_list[i].end - sorted_peaks_list[i+1].start
                if temp_overlap > cur_max_overlap: 
                    cur_max_overlap = temp_overlap
                    merge_pter = i+1

            # then, if the max_overlap you find is >= specified max_overlap then merge. true for first peak always once inputted, but maybe not true for peaks post-merging 
            if cur_max_overlap >= max_overlap:

                # new midpt is: weighted avg of midpts of the peaks. get new params for the new peak
                cur_interval, interval_to_merge = (sorted_peaks_list[merge_pter-1].start, sorted_peaks_list[merge_pter-1].end), (sorted_peaks_list[merge_pter].start, sorted_peaks_list[merge_pter].end)
                cur_mid, next_mid = self.calculate_mid(cur_interval[0], cur_interval[1]), self.calculate_mid(interval_to_merge[0], interval_to_merge[1])
                a, b = sorted_peaks_list[merge_pter-1].count, sorted_peaks_list[merge_pter].count
                mid = (a * cur_mid + b * next_mid) / (a + b)

                # update the new merged peak 
                sorted_peaks_list[merge_pter-1].update_counts(a, b)
                sorted_peaks_list[merge_pter-1].extend_coords(mid, sequence_length)

                # ~ iterate by popping off list 
                peaks_to_merge.pop()
            
            # ~ keep iterating. Should modify sorted_intervals in place 

# merge the peaks, in each chromosome, so that they are 600 bp long. followed procedure from the paper 
def greedily_merge_peaks(max_overlap, file_dict, row_dict, sequence_length):
    
    # 0. for each chromosome => 
    # 1. iterate thru peaks, extending each one to +/- y bp from midpoint. 
    # 2. if there are overlaps of bp >= x, merge them, otherwise leave the newly extended peaks alone 

    for chromosome, all_rows in file_dict.items():

        # sort all intervals by start coordinate, makes the below algorithm more efficient than the while loops make it out to be
        sorted_intervals = sorted(all_rows.keys(), key=lambda x : x[0])
        # create a peaks list for all the intervals in the sorted intervals list 
        sorted_peaks_list = [Peak(interval[0], interval[1], all_rows[(interval[0], interval[1])]) for interval in sorted_intervals]
        # iterate over all the intervals until you encounter a peak that overlaps > 200 bp with the peak you're currently at. the objective is to merge that into our current peak
        peaks_to_merge = []
        
        for i, cur_peak in enumerate(sorted_peaks_list):

            # for every peak just extend it. update entry in sorted_peaks_list 
            cur_peak.extend_coords(cur_peak.mid, sequence_length) 

            # to start off, always have the precaution of having to merge your first peak. means after i = 0, peaks_to_merge is always >1 in length 
            if len(peaks_to_merge) == 0:
                largest_end = cur_peak.end
                peaks_to_merge = [cur_peak]
            else: 
                # update the largest end to whatever peak was previous 
                largest_end = sorted_peaks_list[i-1].end 
                # keep adding peaks to the list of peaks to merge UNTIL you encounter an overlap that is <= the max amount. once u encounter, merge all peaks that overlap >max amt of bp up until then
                if largest_end - max_overlap <= cur_peak.start:
                    cur_peak.merge(peaks_to_merge, max_overlap, sorted_peaks_list, sequence_length)
                else: 
                    peaks_to_merge.append(cur_peak)
        
        # TESTING ~~~
        var = [(peak.start, peak.end) for peak in sorted_peaks_list if peak.end - peak.start != 600]
        if len(var) > 0: 
            print(var)
        
        # update the row_dict w/ the intervals that are newly merged 
        row_dict[chromosome] = sorted_peaks_list

    return row_dict

# notes how active each individual peak, based on its presence across all cell types. for me, age points? 
def input_active_peaks(file_dict:dict, unique_intervals:dict):
    
    row_dict = {}

    # if we're testing by pure equality alone, then just use the existing dict
    for chrom, rows_per_chromosome in file_dict.items(): 
        # should just be 1 interval per assuming non-overlapping windows 
        for interval in rows_per_chromosome.keys():
            rows_per_chromosome[interval] = unique_intervals.get(interval, 1)

    # now populate a new dict that will be the file dict from now on.  
    for chrom, rows_per_chromosome in file_dict.items(): 
        row_dict[chrom] = []
        for interval in rows_per_chromosome.keys(): 
            # chrom => [(start, end), int count]
            row_dict[chrom] += [(interval[0], interval[1])]
    return file_dict, row_dict

if os.path.exists(directory) and os.path.isdir(directory):

    # iterate thru all files 
    unique_intervals = {}
    all_files = []
    for file in os.listdir(directory):
        filepath = os.path.join(directory, file) 
        if os.path.isfile(filepath):
            # generate dictionary for the file 
            file_dict = generate_metrics({}, pb.BedTool(filepath), unique_intervals)
            all_files.append(file_dict)

    # iterate thru again to update the unique interval counts, greedily merge peaks 
    for j in range(len(all_files)): 
        file_dict, row_dict = input_active_peaks(all_files[j], unique_intervals)
        print("FILE:", j+1)
        row_dict = greedily_merge_peaks(max_overlap, file_dict, row_dict, sequence_length)
else: 
    print('Inputted directory filepath does not exist')

# ---------- IGNORE ---------- # 
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