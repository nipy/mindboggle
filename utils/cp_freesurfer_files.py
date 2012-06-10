#!/usr/bin/python

"""
Copy select files from one FreeSurfer subject directory to another.

Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os

atlas_list_file = '/projects/mindboggle/mindboggle/data/atlases/subjects.txt'
atlases_orig_path = '/projects/mindboggle/data/ManualSurfandVolLabels/subjects'
atlases_path = '/projects/mindboggle/mindboggle/data/atlases/subjects_freesurfer'

source_files = os.path.join('surf', '*.orig')
target_dir = 'surf'

# List of atlases, new and original names
f = open(atlas_list_file)
atlases_lines = f.readlines()
atlases = [a.strip("\n").split("\t")[0] for a in atlases_lines if a.strip("\n").split("\t")]
atlases_orig = [a.strip("\n").split("\t")[1] for a in atlases_lines if a.strip("\n").split("\t")]

for i, atlas in enumerate(atlases):
    atlas_orig = atlases_orig[i]

    what_to_copy = os.path.join(atlases_orig_path, atlas_orig, source_files)
    where_to_copy = os.path.join(atlases_path, atlas, target_dir)

    c = 'cp ' + what_to_copy + ' ' + where_to_copy
    print(c)
    os.system(c)

