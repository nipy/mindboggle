# Use FreeSurfer's tksurfer to visualize .annot surface mesh data

import os

subject = 'MMRR-21-1'
hemisphere = 'lh'
surface = 'pial'
annotname = 'labels.max'
colortable = '/projects/mindboggle/mindboggle/data/atlases/atlas_color_LUT.txt'

args = ['tksurfer', subject, hemisphere, surface, '-annotation', annotname,
        '-colortable', colortable]
c = ' '.join(args)
print(c)
os.system(c)
