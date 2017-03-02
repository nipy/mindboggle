#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

colormap_filename = 'colors.npy'
labels_filename = 'IDs.npy'
lut_filename = 'colors.lut'

labels = np.load(labels_filename)
colors = np.load(colormap_filename)
colors = colors * 255

contents = list()
contents.append('<LUT>')
contents.append('256\t# Size of LUT Arrays')

for i in range(256):
    if i in labels:
        idx = np.where(labels == i)
        c = colors[idx].tolist()[0]
    else:
        c = [0.0, 0.0, 0.0]
    c_str = [str(cc) for cc in c]
    line = '\t'.join([str(i), '1.0', *c_str])
    contents.append(line)

with open(lut_filename, 'w') as file:
    file.write('\n'.join(contents))
