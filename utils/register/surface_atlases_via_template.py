#!/usr/bin/python
"""
Use FreeSurfer to register multiple atlas subjects' brain surface labels
to a target subject's brain surface via an average template.

Command: python <this file name>
                <FreeSurfer subjects directory>
                <FreeSurfer subject name> <output directory>

Example: python surface_labels_via_template.py
                /Applications/freesurfer/subjects/ bert output/ 

This program uses FreeSurfer's mris_surf2surf, which
resamples one cortical surface onto another, and assumes that the atlas
directories have been placed in the FreeSurfer subjects directory.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import sys
import os

# Check inputs
if len(sys.argv) == 4:
    subjects_path = sys.argv[1] + '/'
    subject_name = sys.argv[2]
    output_path = sys.argv[3] + '/'
else:
    print("Please check your command-line arguments.")
    sys.exit()

atlases_path = subjects_path
atlas_list_file = output_path + 'atlas_list.txt'
annot_file = 'aparcNMMjt.annot'

subject_path = subjects_path + subject_name
if os.path.exists(subject_path) and\
   os.path.exists(atlases_path) and \
   os.path.exists(atlas_list_file) and \
   os.path.exists(output_path):
    pass
else:
    print(subject_path + ', ' + atlases_path + ', ' + atlas_list_file + \
          ', ' + output_path + " doesn't exist.")
    sys.exit()

# Get list of atlas subjects from a file
f = open(atlas_list_file)
atlas_list = f.readlines()
for atlas_line in atlas_list:
    atlas_name = atlas_line.strip("\n")
    atlas_path = atlases_path + atlas_name
    if os.path.exists(atlas_path):
        pass
    else:
        print(atlas_path + " doesn't exist.")
        sys.exit()

    # For each hemisphere
    for h in ['lh']:
        registration_name = '.sphere_to_template.reg'
        output_annot = output_path + h + '.' + atlas_name + '_to_' + \
                       subject_name + '_' + annot_file

        atlas_to_template_reg = atlas_path + '/surf/' + h + '.sphere.reg'
        subject_to_template_reg = subject_path + '/surf/' + h + '.sphere.reg'
        if os.path.exists(atlas_to_template_reg) and \
           os.path.exists(subject_to_template_reg):
            pass
        else:
            print(atlas_to_template_reg + ' or ' + subject_to_template_reg + \
                  " doesn't exist.")
            sys.exit()

        # Transform atlas labels to the target
        args = ['mri_surf2surf',
                '--hemi', h,
                '--srcsubject', atlas_name,
                '--trgsubject', subject_name,
                '--sval-annot', annot_file,
                '--tval', output_annot,
                '--srcsurfreg', registration_name,
                '--trgsurfreg', registration_name]
        print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);
