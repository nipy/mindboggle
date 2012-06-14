#!/usr/bin/python
"""
Use ANTs to transform an atlas brain image volume's labels
to a subject's brain image volume via an average template.

OR

Use ANTs to affine transform an atlas brain image volume
to a subject's brain image volume via an average template
and compute image similarity between each atlas and the subject.

Command: python <this file name>
                <"1" for label transforms, "0" for image similarity measures>
                <path to atlases directory>
                <atlas name (directory within atlases directory)>
                <path to subjects directory>
                <subject name (directory within subjects directory)>
                <output path>

Example: python volume_atlases_via_template.py  1  ./atlases/  atlas1
                ./subjects/  subject1  ./output/

This program uses ANTs SyN registration method.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import sys

template = '../templates/volume_templates/KKI_template.nii.gz'

# Check inputs
if len(sys.argv) == 7:
    transform_labels = sys.argv[1] + '/'
    atlases_path = sys.argv[2] + '/'
    atlas_name = sys.argv[3]
    subjects_path = sys.argv[4] + '/'
    subject_name = sys.argv[5]
    output_path = sys.argv[6] + '/'
else:
    print("Please check your command-line arguments.")
    sys.exit()

atlas_path = atlases_path + atlas_name + '/'
subject_path = subjects_path + subject_name + '/'
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

# Nonlinearly transform atlas labels to the target
if transform_labels:
    atlas = atlas_path + 'labels.nii.gz'
    output = output_path + atlas_name + '_to_subject.nii.gz'
    args = ['WarpImageMultiTransform 3', atlas, output, '-R', subject,
            '-i', subject_to_template_reg + 'Affine.txt',
            subject_to_template_reg + 'InverseWarp.nii.gz',
            atlas_to_template_reg + 'Warp.nii.gz',
            atlas_to_template_reg + 'Affine.txt --use-NN']

    cmd = 'date'; os.system(cmd)
    print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);
    cmd = 'date'; os.system(cmd)

# Affine transform atlas image to the target
else:
    atlas = atlas_path + 'image.nii.gz'
    output = output_path + atlas_name + '_to_subject_affine.nii.gz'
    args = ['WarpImageMultiTransform 3', atlas, output, '-R', subject,
            '-i', subject_to_template_reg + 'Affine.txt',
            atlas_to_template_reg + 'Affine.txt']
    print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);

    args = ['ImageSimilarity 3', output, subject]
    print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);
