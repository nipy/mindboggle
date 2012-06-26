#!/usr/bin/python

"""
Use FreeSurfer's tksurfer to visualize .annot surface mesh data

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0
"""

import os, sys

debug = 0

if len(sys.argv) < 6:
    if debug:
        subject = 'MMRR-21-1'
        hemisphere = 'lh'
        surface = 'pial'
        annotname = 'labels.max'
        colortable = '/projects/mindboggle/mindboggle/data/atlases/labels.DKT26.txt'
    else:
        sys.exit('Usage: %s subject hemisphere surface-type annotname colortable' %sys.argv[0])
else:
    subject = sys.argv[1]  #'MMRR-21-1'
    hemisphere = sys.argv[2]  #'lh'
    surface = sys.argv[3]  #'pial'
    annotname = sys.argv[4]  #'labels.max'
    colortable = sys.argv[5]  #'/projects/mindboggle/mindboggle/data/atlases/labels.DKT26.txt'

args = ['tksurfer', subject, hemisphere, surface, '-annotation', annotname,
        '-colortable', colortable]
c = ' '.join(args)
print(c); os.system(c)
