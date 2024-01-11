import numpy as np
import sys
from matplotlib import pyplot as plt 
import pandas as pd

# NOTE: i should probably get rid of the header in the file

# filepath = sys.argv[1]
filepath = 'final_counts_draft.csv'
xy = np.loadtxt(filepath, delimiter=",", dtype='str')
print(np.shape(xy)) # (13894, 5)

setr1_counts = np.array(xy[:,2], dtype='int')
setr2_counts = np.array(xy[:,3], dtype='int')
setr3_counts = np.array(xy[:,4], dtype='int')

# now plot the above counts => distribution of read counts for all conditions 
f, a = plt.subplots(1, 3, figsize=(10,7.5), layout='constrained')
f.suptitle('Read Count Distributions')
print('-------------- setr1 summary stats ------------')
print('max: ', np.max(setr1_counts), 'min: ', np.min(setr1_counts), 'mean: ', np.mean(setr1_counts), 'median', np.median(setr1_counts))
a[0].set_xlabel('Frequency')
a[0].set_ylabel('Read Counts')
a[0].set_title('Setr1 Read Counts')
a[0].hist(setr1_counts, bins=np.linspace(1,100))

print('-------------- setr2 summary stats ------------')
print('max: ', np.max(setr2_counts), 'min: ', np.min(setr2_counts), 'mean: ', np.mean(setr2_counts), 'median', np.median(setr2_counts))
a[1].set_xlabel('Frequency')
a[1].set_ylabel('Read Counts')
a[1].set_title('Setr2 Read Counts')
a[1].hist(setr2_counts, bins=np.linspace(1,100))

print('-------------- setr3 summary stats ------------')
print('max: ', np.max(setr3_counts), 'min: ', np.min(setr3_counts), 'mean: ', np.mean(setr3_counts), 'median', np.median(setr3_counts))
a[2].set_xlabel('Frequency')
a[2].set_ylabel('Read Counts')
a[2].set_title('Setr3 Read Counts')
a[2].hist(setr3_counts, bins=np.linspace(1,100))
plt.show()

