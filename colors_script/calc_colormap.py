#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append('../mindboggle/mio')
import colors
import numpy as np
import os
import pickle

filename = 'label.nii.gz'
matrix_filename = 'matrix.pickle'
colormap_filename = 'colormap.npy'

if not os.path.isfile(matrix_filename):
    print('finding adjacency maps')
    IDs, matrix, output = colors.label_adjacency_matrix(filename, ignore_values=[])
    pickle.dump(matrix, open(matrix_filename, 'wb'))
    np.save('IDs', IDs)
else:
    IDs = np.load('IDs.npy')
    matrix = pickle.load(open(matrix_filename, 'rb'))

matrix = matrix.as_matrix()
matrix = matrix[:, 1:]

print('finding colormap')
if not os.path.isfile(colormap_filename):
    colormap = colors.distinguishable_colors(ncolors=len(IDs),
                                            backgrounds=[[0,0,0],[1,1,1]],
                                            save_csv=False)
    np.save(colormap_filename, colormap)
else:
    colormap = np.load(colormap_filename)

print(colormap.shape)

colors = colors.group_colors(colormap, 'cerebellum', IDs=IDs,
                             adjacency_matrix=matrix)
np.save('colors', colors)
