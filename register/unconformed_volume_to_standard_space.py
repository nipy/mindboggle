#!/usr/bin/python

"""
Apply FreeSurfer's "Talairach" transform matrix to take an unconformed volume to a standard space.

Command: python <this file name> <input nifti volume> <path to FreeSurfer subject> 
                <output directory> <output name stem> <interpolation>

Example: python unconformed_volume_to_standard_space.py
                subject1.nii.gz /Applications/freesurfer/subjects/bert output/ bert nearest

See http://surfer.nmr.mgh.harvard.edu/fswiki/mri_convert

Batch example: Run transform_unconformed_volume_to_standard_space on a directory:
  import os
  output_path = 'transformed_brainvisa_fundus_volumes_Perrot62'
  subjects_path = '/home/arno/Data/Brains/Perrot62_sulci/freesurfer5.1_output_plus_surface_features/'
  volumes_path = '/home/arno/Data/Brains/Perrot62_sulci/manually_labeled_brainvisa_fundi/sulci_volumes/'
  subjects = os.listdir(subjects_path)
  volumes = os.listdir(volumes_path)
  for i,volume in enumerate(volumes):
    args = ['python transform_unconformed_volume_to_standard_space.py',
            volumes_path+volume, subjects_path+subjects[i], output_path, volume, 'nearest']
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys

# Check inputs
if len(sys.argv) < 6:
    print("Usage: python transform_unconformed_volume_to_standard_space.py <input nifti volume> <path to FreeSurfer subject> <output directory> <output name stem>")
    exit(-1)
else:
    native_volume_nii = sys.argv[1]
    subject_path = sys.argv[2] + '/'
    output_path = sys.argv[3] + '/'
    output_name = sys.argv[4]
    interpolation = sys.argv[5]  # Ex: nearest, interpolate, sinc, cubic

# Apply FreeSurfer's "Talairach" transform matrix to take a conformed volume to standard space:
if os.path.exists(native_volume_nii):
    pass
else:
    print(native_volume_nii + " doesn't exist.")
    exit(-1)

conformed_volume_mgz = output_path + output_name + '.conformed.mgz'
transformed_volume_nii = output_path + output_name + '.transformed.nii.gz'
xfm = subject_path + 'mri/transforms/talairach.xfm'

args = ['mri_convert --conform -rt', interpolation, native_volume_nii, conformed_volume_mgz]
print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
args = ['mri_convert --apply_transform', xfm, '-rt',interpolation, conformed_volume_mgz, transformed_volume_nii]
print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
