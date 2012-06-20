#!/usr/bin/python

"""
Copy select files from one FreeSurfer subject directory to another.

Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import re

atlas_list_file = '/projects/mindboggle/mindboggle/data/atlases/atlas_list_new_old.txt'
atlases_orig_path = '/home/arno/Data/Brains/Mindboggle100_Neuromorphometrics/subjects'
atlases_path = '/projects/mindboggle/mindboggle/data/atlases/atlases_freesurfer'

source_files = os.path.join('mri', 'aseg.mgz')
target_dir = 'mri'

# List of atlases, new and original names
f = open(atlas_list_file)
atlases_lines = f.readlines()

atlases = []
atlases_orig = []
for atlas_line in atlases_lines:
    atlases.append(re.findall(r'\S+', atlas_line)[0])
    atlases_orig.append(re.findall(r'\S+', atlas_line)[1])

for i, atlas in enumerate(atlases):
    atlas_orig = atlases_orig[i]

    what_to_copy = os.path.join(atlases_orig_path, atlas_orig, source_files)
    where_to_copy = os.path.join(atlases_path, atlas, target_dir)

    c = 'cp ' + what_to_copy + ' ' + where_to_copy
    print(c)
    os.system(c)

