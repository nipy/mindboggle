#!/usr/bin/env python
"""
Visualize .label labeled surface mesh data.

Comment out the top line (below the header) of the .label file
containing the number of rows, to load an appropriately sized numpy array.

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
"""

import sys
import numpy as np
from matplotlib import pyplot as plt

print('NOTE: ')
print('Comment out the top line (below the header) of the .label file')
print('containing the number of rows, to load an appropriately sized numpy array.')

if len(sys.argv) < 2:
    sys.exit('Usage: %s label_file' %sys.argv[0])
else:
    label_file = sys.argv[1]

a = np.loadtxt(label_file)

if a.ndim == 2:

    x = a[:,1]
    y = a[:,2]
    z = a[:,3]

    fig = plt.figure()
    ax = fig.add_subplot(111) #, projection='3d')
    ax.plot(x,y,z, 'k.')
    plt.show()

else:
    sys.exit("The label array is the wrong shape." +\
    "Remember to comment out the top line (below the header) of the .label file" +\
    "containing the number of rows.")
