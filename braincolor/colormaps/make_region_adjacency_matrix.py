#! /usr/bin/env python
"""                                                      
http://www.braincolor.org

Create a region adjacency matrix.

Input:
1. Labeled image volume, where each region has an index value.
2. Text file, with each line containing:
region abbreviation, region name, and group number 
(group number 1...N indicates inclusion within a color family): 
----------------------------------------
3 8 Cb "L Cerebellum Exterior"
4 8 Cb "R Cerebellum Exterior"
15 7 Hip "L Hippocampus"
16 7 Hip "R Hippocampus"
...
----------------------------------------

Output:
Columns 1-3 will contain the same information as the input file above.
The remaining columns and rows will contain the automatically generated adjacency matrix: 
----------------------------------------
ID Abbr Name Group Cb Cb Hip Hip ...
3 8 Cb "L Cerebellum Exterior" 0 0 0 ...
4 8 Cb "R Cerebellum Exterior" 0 0 0 ...
15 7 Hip "L Hippocampus" 0 0 0 ...
16 7 Hip "R Hippocampus" 0 0 0 ...
...
----------------------------------------

(c) Copyright 2011 . arno klein . arno@binarybottle.com . MIT license
"""

from subprocess import call
from nibabel import load
from numpy import nonzero, ravel, take, size, unique, zeros, intersect1d
#from scipy import ndimage

# Region adjacency table
in_dir = 'input/CUMC12/'
in_volume = in_dir + 'm1.nii.gz'
in_table = in_dir + 'CUMC12_labels_regions.txt'
out_temp = in_dir + 'temp.nii.gz'
out_table = in_dir + 'region_adjacency_matrix.csv'
col_ID = 0  # column with region abbreviations (not used by the program)
col_group = 1  # column with region group numbers
col_abbr = 2  # column with region abbreviations
col_name = 3  # column with full region names (not used by the program)
col_start_data = 4  # first column with data
row_start_data = 0  # first row with data
everyother = 2  # use <everyother> alternate row(s);
                # set to 2 for redundant labels across, e.g. brain hemispheres
                
f = open(in_table,'r')
table = f.readlines()
f.close()
f = open(out_table,'w')
#f.close()
#f = open(out_table,'a')
L = load(in_volume).get_data()

labels = []
for row in table:
    labels.append(int(row.split(" ")[col_ID]))
max_label = max(labels)

for row in table:
    ID = row.split(" ")[col_ID]

    ## If scipy is installed, don't use Convert3d:
    # B = L * (ndimage.binary_dilation(L==int(ID)) - (L==int(ID)))
    ## If scipy is not installed, use Convert3d:
    args = " ".join(['c3d',in_volume,'-threshold',ID,ID,'1 0 -o', out_temp])
    print(args);  p = call(args, shell="True")  
    args = " ".join(['c3d',out_temp,'-dilate 1 1x1x1mm -o', out_temp])
    print(args);  p = call(args, shell="True")  
    T = load(out_temp).get_data()
    B = L * (T - (L==int(ID)))

    neighbor_labels = unique(ravel(B))
    neighbor_labels = [int(s) for s in neighbor_labels if s>0 and s in labels and s%2==(int(ID)%2)]
    adjacency_vector = zeros(max_label)
    for ivec in neighbor_labels:
        adjacency_vector[ivec-1] = 1
    adjacency_vector = take(adjacency_vector,[int(s-1) for s in labels])
    adjacency_vector = " ".join([str(s) for s in adjacency_vector])
    if len(adjacency_vector) > 0:
        neighbor_labels = " ".join([str(s) for s in neighbor_labels])
        print(" ".join([row.strip("\n"), neighbor_labels]))
        print(" ".join([row.strip("\n"), adjacency_vector]))
        f.write(" ".join([row.strip("\n"), adjacency_vector, "\n"]))

f.close()
