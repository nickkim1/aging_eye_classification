
# NOTE: installed bedtools w/ brew locally 

# wget command: 
# wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E003-DNase.hotspot.fdr0.01.peaks.v2.bed.gz

# This is only when the window size is SHORTER than the actual window size

# Sort input files with bedtools  
for file in ../0_data/1_exp_data/bed/raw/*.bed; do 
    filename=$(basename -- "$file")
    without_extension=${filename%.bed}
    bedtools sort -i $file > ../0_data/1_exp_data/bed/sorted/$without_extension.sorted.bed
done 

# Partition into non-overlapping windows 
for file in ../0_data/1_exp_data/bed/sorted/*.sorted.bed; do
    filename=$(basename -- "$file")
    without_extension=${filename%.sorted.bed}
    bedtools makewindows -b $file -w 20 > ../0_data/1_exp_data/bed/binned/$without_extension.binned.bed
done 

# Intersect to get proper read counts within each window 
for file in ../0_data/1_exp_data/bed/binned/*.binned.bed; do
    filename=$(basename -- "$file")
    without_extension=${filename%.binned.bed}
    bedtools intersect -a $file -b ../0_data/1_exp_data/bed/raw/$without_extension.bed -c > ../0_data/1_exp_data/bed/intersected/$without_extension.intersected.bed
done 