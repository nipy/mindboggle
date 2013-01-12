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

if len(sys.argv) < 2:
    #sys.exit('Usage: %s label-file-name' %sys.argv[0])
    label_file = '/projects/Mindboggle/output/workspace/Mindboggle/Fill_volume_prep/_hemi_lh_subject_HLN-12-3/Write_label_files/mapflow/_Write_label_files6/lh.inferiortemporal.label'
    #label_file = '/projects/Mindboggle/output/workspace/Mindboggle/Fill_volume_prep/_hemi_lh_subject_HLN-12-3/Write_label_files/mapflow/_Write_label_files30/lh.insula.label'
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
