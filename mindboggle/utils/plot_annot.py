#!/usr/bin/python

"""
Use FreeSurfer's tksurfer to visualize .annot surface mesh data

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
"""

import os, sys

debug = 0

if len(sys.argv) < 6:
    if debug:
        subject = 'HLN-12-3'
        hemisphere = 'lh'
        surface = 'pial'
        annotname = 'labels.DKT31.manual'
        colortable = '/mindboggle/info/labels.surface.DKT31.txt'
    else:
        sys.exit('Usage: %s subject hemisphere surface-type annotname colortable' %sys.argv[0])
else:
    subject = sys.argv[1]
    hemisphere = sys.argv[2]
    surface = sys.argv[3]
    annotname = sys.argv[4]
    colortable = sys.argv[5]

args = ['tksurfer', subject, hemisphere, surface, '-annotation', annotname,
        '-colortable', colortable]
c = ' '.join(args)
print(c); os.system(c)
