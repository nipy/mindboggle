#!/usr/bin/python
"""
Use FreeSurfer's tksurfer to visualize .annot surface mesh data

Input: VTK surface mesh file name

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
"""

import os, sys

if len(sys.argv) < 2:
    sys.exit('Usage: %s vtk-file-name' %sys.argv[0])
else:
    vtk_file = sys.argv[1]

args = ['mayavi2', '-d' , vtk_file, '-m Surface']
c = ' '.join(args)
print(c); os.system(c)
