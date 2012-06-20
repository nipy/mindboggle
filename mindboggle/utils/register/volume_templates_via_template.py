#!/usr/bin/python
"""
Use ANTs to register a subject's multiple atlas subjects' brain volume labels
to a target subject's brain volume via an average template.

Command: python <this file name> <path to subject directory>
                                 <template file> <output path>

Example: python volume_to_template.py
                ./subjects/  bert  ./templates_ants/KKI.nii.gz

This program uses ANTs SyN registration method.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import sys

template = '../templates/volume_templates/KKI_template.nii.gz'

# Check inputs
if len(sys.argv) == 6:
    atlases_path = sys.argv[1] + '/'
    atlas_name = sys.argv[2]
    subjects_path = sys.argv[3] + '/'
    subject_name = sys.argv[4]
    output_path = sys.argv[5] + '/'
else:
    print("Please check your command-line arguments.")
    sys.exit()

atlas_path = atlases_path + atlas_name
subject_path = subjects_path + subject_name
if os.path.exists(atlas_path) and\
   os.path.exists(subject_path) and\
   os.path.exists(output_path):
    pass
else:
    print(atlas_path + ', ' + subject_path + ', or ' + output_path + \
          " doesn't exist.")
    sys.exit()

atlas_to_template_reg = atlas_path + 'registered_to_template_'
subject_to_template_reg = subject_path + 'registered_to_template_'
atlas_reg_test = atlas_to_template_reg + 'Warp.nii.gz'
subject_reg_test = subject_to_template_reg + 'Warp.nii.gz'
if os.path.exists(atlas_reg_test) and os.path.exists(subject_reg_test):
    pass
else:
    print(atlas_reg_test + ' or ' + subject_reg_test + " doesn't exist.")
    sys.exit()

# Transform atlas labels to the target
output = output_path + atlas_name + '_to_subject.nii.gz'
args = ['WarpImageMultiTransform 3', atlas, output, '-R', subject,
        '-i', subject_to_template_reg + 'Affine.txt',
        subject_to_template_reg + 'InverseWarp.nii.gz', 
        atlas_to_template_reg + 'Warp.nii.gz',
        atlas_to_template_reg + 'Affine.txt --use-NN']

cmd = 'date'; os.system(cmd)
print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);
cmd = 'date'; os.system(cmd)
