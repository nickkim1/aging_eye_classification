import pybedtools as pb
import argparse as ap
import os 

# merge peaks for atac-seq based data 
parser = ap.ArgumentParser(prog='This program processes in ATACseq data in BED format. Bins peaks into 600 bp bins')
parser.add_argument('-directory', type=str, help='Input directory holding BED files to process, in comma delimited format', nargs=1)
parser.add_argument('-maxoverlap', type=int, help='Maximum overlap between consecutive intervals in the BED file')
parser.add_argument('-sequencelength', type=int, help="How long you want input sequences to the model to be")
args = vars(parser.parse_args())

# set global input arguments to be used in the script
directory = args['directory'][0]
max_overlap = args['maxoverlap'][0]
sequence_length = args['sequencelength'][0]

# ------------------------------------------------------------------------------- # 
# Methods for generating merged results 
# ------------------------------------------------------------------------------- # 

# generate active input metrics 
def generate_metrics(file_dict:dict, file, unique_intervals:dict):
    for row in file: 
        file_dict[row.chrom] = {} 
        file_dict[row.chrom][(int(row[1]), int(row[2]))] = [int(float(row[4]))]
        unique_intervals[(int(row[1]), int(row[2]))] = 1 + unique_intervals.get((int(row[1]), int(row[2])), 0)
    return file_dict

# merge the peaks, in each chromosome, so that they are 600 bp long. followed procedure from the paper 
def greedily_merge_peaks(max_overlap, file_dict, row_dict, sequence_length):
    
    # 0. for each chromosome => 
    # 1. iterate thru peaks, extending each one to +/- y bp from midpoint. 
    # 2. if there are overlaps of bp >= x, merge them, otherwise leave the newly extended peaks alone 

    for chromosome, all_intervals in row_dict.items():

        # sort all intervals by start coordinate, makes the below algorithm more efficient than the while loops make it out to be
        sorted_intervals = sorted(all_intervals, key=lambda x : x[0])

        # iterate over all the intervals until you encounter a peak that overlaps > 200 bp with the peak you're currently at. do not merge that one
        peaks_to_merge = []
        for i, interval in enumerate(sorted_intervals): 

            # start with a peak, and extend it 
            cur_peak = Peak(interval[0], interval[1])
            cur_peak.extend_coords(sequence_length) 
            largest_end = 0

            # to start off, always have the precaution of having to merge your first peak. means after i = 0, peaks_to_merge is always >1 in length 
            if len(peaks_to_merge) == 0:
                peaks_to_merge = [cur_peak]
            else: 
                # if you've extended the peak, but it overlaps with a previous peak with x amt of bp, then merge it 
                if largest_end - cur_peak.start >= max_overlap: 

                    # start merging peaks previous to this point. exit condition mentioned above
                    cur_peak.merge(cur_peak, peaks_to_merge, max_overlap, i, sorted_intervals, file_dict, sequence_length, chromosome)
                else: 
                    # max overlap isn't exceeded, so just append the peaks to this list. you'll have a separate method that iterates through this list and resets it once you're done with it
                    largest_end = max(largest_end, cur_peak.end)
                    peaks_to_merge.append(cur_peak)
        
        # update the row_dict w/ the intervals that are newly merged 
        row_dict[chromosome] = sorted_intervals

    return row_dict

# notes how active each individual peak, based on its presence across all cell types. for me, age points? 
def input_active_peaks(file_dict:dict, unique_intervals:dict):
    row_dict = {}
    # if we're testing by pure equality alone, then just use the existing dict
    for chrom, rows_per_chromosome in file_dict.items(): 
        # should just be 1 interval per assuming non-overlapping windows 
        for interval in rows_per_chromosome.keys():
            rows_per_chromosome[interval][1] = unique_intervals.get(interval, 1)
    # now populate a new dict that will be the file dict from now on. TODO: optimize this for space 
    for chrom, rows_per_chromosome in file_dict.items(): 
        for interval in rows_per_chromosome.keys(): 
            row_dict[chrom] = [interval + rows_per_chromosome[interval]]
    return file_dict, row_dict

def main(): 
    if os.path.exists(directory) and os.path.isdir(directory):
        # iterate thru all files 
        unique_intervals = {}
        all_files = []
        for file in os.listdir(directory):
            # logs each file into a dictionary  
            file_dict = {}
            # check if item is a file 
            filepath = os.path.join(directory, file)
            if os.path.isfile(filepath):
                # generate dictionary for the file 
                file_dict = generate_metrics(file_dict, pb.BedTool(filepath), unique_intervals)
                all_files.append(file_dict)
        # iterate thru again to update the unique interval counts, greedily merge peaks 
        for j in range(len(all_files)): 
            file_dict, row_dict = input_active_peaks(all_files[j], unique_intervals)
            row_dict = greedily_merge_peaks(max_overlap, file_dict, row_dict, sequence_length)
    else: 
        print('Inputted directory filepath does not exist')


# ------------------------------------------------------------------------------- # 
# Class definition for Peak data 
# ------------------------------------------------------------------------------- # 
class Peak():
    def __init__(self, start, end): 
       self.start = start
       self.end = end
       self.mid = self.calculate_mid(self.start, self,end)
    
    def extend_coords(self, mid, sequence_length): 
        self.start = mid - sequence_length / 2
        self.end = mid + sequence_length / 2
        # TODO: have some exception handling for the end of the chromosome 

    def calculate_mid(self, s, e):
        return (s + e) / 2

    # function that merges any given new peak
    def merge(self, cur_peak, peaks_to_merge:list, max_overlap, cur_pter, sorted_intervals, file_dict:dict, sequence_length, chromosome): 

        # take into account your current peak, merge based on midpoints
        cur_max_overlap = max_overlap

        # the condition before guaranteed overlap <= 200 bp using UNMERGED peaks. 
        # but after you merge peaks, you will likely have additional overlap because you're extending each peak +/- 300 bp
        
        while len(peaks_to_merge) > 1 and cur_max_overlap >= max_overlap: 
            
            # first, try to find max_overlap in peaks to merge. in the beginning, the max_overlap will be guaranteed <= cur_max_overlap. but as you keep merging that might not be the case and you might not need to merge anymore
            next_pter = 0
            for i in range(len(peaks_to_merge)): 
                # if negative, that means NO overlap. if positive, then YES overlap. hence greater than sign 
                temp_overlap = peaks_to_merge[i].end - peaks_to_merge[i+1].start 
                if temp_overlap > cur_max_overlap: 
                    cur_max_overlap = temp_overlap
                    next_pter = i
            
            # then, if the max_overlap you find is >= specified max_overlap then merge. true for first peak always once inputted, but maybe not true for peaks post-merging 
            if cur_max_overlap >= max_overlap: 

                # new midpt is? weighted avg of midpts of the peaks. update new peak with this information 
                cur_interval, next_interval = (cur_peak.start, cur_peak.end), (peaks_to_merge[next_pter].start, peaks_to_merge[next_pter].end)
                cur_mid, next_mid = self.calculate_mid(cur_interval[0], cur_interval[1]), self.calculate_mid(next_interval[0], next_interval[1])
                a, b = file_dict[chromosome][cur_interval][1], file_dict[chromosome][next_interval][1]
                mid = (a * cur_mid + b * next_mid) / (a + b)
                start, end = self.extend_coords(mid, sequence_length)

                # update the entry in sorted_intervals 
                sorted_intervals[cur_pter] = [start, end, a + b]

                # remove the next entry from peaks_to_merge once you merge it into the current peak
                peaks_to_merge.remove(peaks_to_merge[next_pter])
            
            # ~ keep iterating. Should modify sorted_intervals in place 

# run the mainline 
if __name__ is 'main':
    main()