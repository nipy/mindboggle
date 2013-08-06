#!/usr/bin/env python
"""
Use FreeSurfer's tksurfer to visualize .annot surface mesh data

    Example arguments
    -----------------
    subject = 'Twins-2-1'
    hemisphere = 'lh'
    surface = 'pial'
    annotname = 'labels.DKT25.manual'
    colortable = '../info/labels.surface.DKT31.txt'

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License
"""
import os, sys
from mindboggle.utils.utils import execute

if len(sys.argv) < 6:
    sys.exit('Usage: %s subject hemisphere surface-type annotname colortable' %sys.argv[0])
else:
    subject = sys.argv[1]
    hemisphere = sys.argv[2]
    surface = sys.argv[3]
    annotname = sys.argv[4]
    colortable = sys.argv[5]

cmd = ['tksurfer', subject, hemisphere, surface, '-annotation', annotname,
       '-colortable', colortable]
execute(cmd)