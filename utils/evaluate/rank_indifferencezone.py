"""
 Indifference-zone ranking algorithm
 
 This program generates indifference-zone ranking
 to find the top-ranking method (column).
 The notion of an "indifference zone" determining
 a practically significant difference between sample means
 comes from Bechhofer (1954). 

 python rank_indifferencezone.py <file> <int> <string1> <string2>...
 argument 1: table file name
 argument 2: set to 1 to plot notched box plots (optional)
 argument 3: method names (optional)

 EXAMPLE:
 python rank_indifferencezone.py results/table.txt 1 method1 method2 method3...

 INPUT: Text table (N rows x (M+1) columns).
 Column 1 contains a delta for each row, where delta is equal to
 the smallest practically significant difference for that row's data.
 Each of the remaining M columns contains measures/samples for a given method.

 OUTPUT: Mean ranks (1-D array of length M), standard deviations.
 Interpretation: the method with the highest mean rank has performed
 practically significantly better than the other methods if its mean rank
 is at least a standard deviation higher than the others.

 ------------------------------------------------------------------------------
 Copyright (c) 2010 . arno klein . arno@mindboggle.info

 MIT License:
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 ------------------------------------------------------------------------------
"""

import sys
import numpy as np

# Inputs:
input_table = sys.argv[1]
average_zeros = 1 # 1: include trivial (zero) vectors when averaging

if len(sys.argv) > 2:
    import pylab as plt
    plot_data = sys.argv[2]
else:
    plot_data = 0

if len(sys.argv) > 3:
    print_methods = 1
    methods = []
    for iarg in range(3,len(sys.argv)):
        methods.append(sys.argv[iarg])
else:
    print_methods = 0

# Load table:
print(" ")
print("Table:")
print(input_table)
print(" ")
table = np.loadtxt(input_table)
#table = np.reshape(table,[760,3])

# Initialize:
deltas    = table[:,0]
data      = table[:,1::]
nmeasures = np.shape(data)[0]
nmethods  = np.shape(data)[1]
Ranks     = np.zeros((nmethods,nmethods-1,nmeasures))
imethods  = range(nmethods)

# Loop through rows of data to rank methods:
for irow in range(nmeasures):

    row = data[irow,:]
    delta = deltas[irow]
    #delta = 0.1
    # Subtract delta from each row and set matrix R elements to -1, 0, or 1
    # when a value for a method is at least delta less than, within delta of,
    # or at least delta greater than the values for the other methods.
    R = np.zeros((nmethods,nmethods-1))
    for i in range(nmethods):
        truncated_row = row[[j for j in imethods if j!=i]]
        R[i,(truncated_row - delta) > row[i]] = -1
        R[i,(truncated_row + delta) < row[i]] =  1

    Ranks[:,:,irow] = R

# Compute mean ranks and mean ranks within 1 standard deviation of the top:
mean_data  = np.mean(data,axis=0)
std_data   = np.std(data,axis=0,ddof=1)
mean_ranks = np.mean(np.mean(Ranks,2),1)
if Ranks.shape[1] > 1:
    std_ranks = np.std(np.mean(Ranks,2),1,ddof=1)
else:
    std_ranks = np.std(Ranks,2,ddof=1)
top_rank_members = mean_ranks + std_ranks[mean_ranks.argmax()] > np.max(mean_ranks)
top_rank_members = [k+1 for k in imethods if top_rank_members[k]==True and k!=mean_ranks.argmax()]

# Print results:
if print_methods:
    print("Methods:")
    print(methods)
    print(" ")
print("Method (column) means:")
print(mean_data)
print("Method (column) SDs:")
print(std_data)

print("Rank means:")
print(mean_ranks)
print("Rank SDs:")
print(std_ranks)
print(" ")

if print_methods:
    print(methods[mean_ranks.argmax()] + " obtained the maximum mean rank.\n")
else:
    print("Method (column) " + str(mean_ranks.argmax()+1) + " obtained the maximum value.\n")

if any(top_rank_members):
    print("Methods within one standard deviation of the top rank:  " + str(top_rank_members))
else:
    print("No methods were within one standard deviation of the top rank.")
print("")


# Plot results:
if plot_data:
    plt.figure()
    plt.boxplot(data,1)
    plt.title('Measures per method')
    plt.xlabel('Methods')
    plt.ylabel('Measures')
    if print_methods:
        plt.xticks(np.arange(1,nmethods+1), methods, rotation=-15)
