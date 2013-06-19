#!/usr/bin/env python
"""
Post-process Mindboggle-101 volume images for distribution,
using Mindboggle, FreeSurfer, and FSL tools:

    # Convert label volume from FreeSurfer to original space
    # Extract brain by masking with manual cortical
      and automated subcortical labels
    # Remove subcortical labels
    # Convert DKT31 to DKT25 labels
    # Affine register T1-weighted brain to MNI152 brain
    # Transfer whole-head images with affine transform
    # Transfer labeled images with affine transform
      (and nearest neighbor interpolation)

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2013  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os

# Paths, template, and label conversion files
mb101_path = '/hd2/Lab/Brains/Mindboggle101/'
mb_info_path = '/projects/Mindboggle/mindboggle/mindboggle/info/'
template = mb101_path+'MNI152space/MNI152_T1_1mm_brain.nii.gz'
relabel_file = os.path.join(mb_info_path, 'labels.volume.DKT31to25.txt')

# Loop through subjects
list_file = mb_info_path + 'atlases101.txt'
fid = open(list_file, 'r')
subjects = fid.readlines()
subjects = [''.join(x.split()) for x in subjects]
for subject in subjects:

    print(">>> Process subject: {0}...".format(subject))
    subject_path = mb101_path + 'subjects/' + subject + '/mri/'

    # Identify original files
    full_labels_orig = subject_path+'aparcNMMjt+aseg.nii.gz'
    head = subject_path+'t1weighted.nii.gz'

    # Name all output files
    local_labels0 = 'labels.DKT31.manual+aseg.nii.gz'
    full_labels = subject_path + local_labels0
    local_labels = 'labels.DKT31.manual.nii.gz'
    DKT31_labels = subject_path + local_labels
    DKT25_labels = subject_path + 'labels.DKT25.manual.nii.gz'
    brain = subject_path+'t1weighted_brain.nii.gz'
    xfm_matrix = subject_path+'t1weighted_brain.MNI152.mat'
    xfm_brain = subject_path+'t1weighted_brain.MNI152.nii.gz'
    xfm_head = subject_path+'t1weighted.MNI152.nii.gz'
    xfm_DKT25 = subject_path+'labels.DKT25.manual.MNI152.nii.gz'
    xfm_DKT31 = subject_path+'labels.DKT31.manual.MNI152.nii.gz'
    xfm_DKT31aseg = subject_path+'labels.DKT31.manual+aseg.MNI152.nii.gz'

    # Convert label volume from FreeSurfer 'unconformed' to original space
    print("Convert label volume from FreeSurfer 'unconformed' to original space...")
    #if 'OASIS-TRT-20-' in subject or 'NKI-TRT-20-' in subject:
    cmd = ' '.join(['mri_vol2vol --mov', full_labels_orig, '--targ', head,
                    '--regheader --o', full_labels])
    #cmd = ' '.join(['mri_convert -rl', head, '-rt nearest',
    #                full_labels_orig, full_labels])
    print(cmd); os.system(cmd)

    # Extract brain by masking with labels using FreeSurfer
    print("Extract brain by masking with labels using FreeSurfer...")
    cmd = ' '.join(['mri_vol2vol --mov', full_labels, '--targ', head,
                    '--o temp.nii.gz --regheader'])
    print(cmd); os.system(cmd)
    cmd = ' '.join(['mri_mask', head, 'temp.nii.gz', brain])
    #cmd = ' '.join(['/usr/bin/fsl4.1-fslmaths', head, '-mas', full_labels, brain])
    print(cmd); os.system(cmd)

    # Remove subcortical labels
    print("Remove subcortical labels...")
    from mindboggle.label.relabel import remove_volume_labels
    labels_to_remove = range(1,300) # Remove noncortical (+aseg) labels
    labels_to_remove.extend([1000,1001,2000,2001])
    remove_volume_labels(full_labels, labels_to_remove)
    cmd = ' '.join(['mv', local_labels0, DKT31_labels])
    print(cmd); os.system(cmd)

    # Convert DKT31 to DKT25 labels
    print("Convert DKT31 to DKT25 labels...")
    from mindboggle.utils.io_file import read_columns
    from mindboggle.label.relabel import relabel_volume
    old_labels, new_labels = read_columns(relabel_file, 2)
    relabel_volume(DKT31_labels, old_labels, new_labels)
    cmd = ' '.join(['mv', local_labels, DKT25_labels])
    print(cmd); os.system(cmd)

    # Affine register T1-weighted brain to MNI152 brain using FSL's flirt
    print("Affine register T1-weighted brain to MNI152 brain using FSL's flirt...")
    cmd = ' '.join(['flirt', '-in', brain, '-ref', template,
                    '-out', xfm_brain, '-omat', xfm_matrix])
    print(cmd); os.system(cmd)

    # Transfer whole-head images with affine transform using FSL's flirt
    print("Transfer whole-head images with affine transform using FSL's flirt...")
    cmd = ' '.join(['flirt', '-in', head, '-ref', template,
                    '-applyxfm -init', xfm_matrix, '-out', xfm_head])
    print(cmd); os.system(cmd)

    # Transfer labeled images with affine transform (and nearest neighbor interpolation)
    print("Transfer labeled images with affine transform "
          "(and nearest neighbor interpolation)...")
    cmd = ' '.join(['flirt', '-in', DKT25_labels, '-ref', template,
                    '-applyxfm -init', xfm_matrix,
                    '-interp nearestneighbour -out', xfm_DKT25])
    print(cmd); os.system(cmd)

    cmd = ' '.join(['flirt', '-in', DKT31_labels, '-ref', template,
                    '-applyxfm -init', xfm_matrix,
                    '-interp nearestneighbour -out', xfm_DKT31])
    print(cmd); os.system(cmd)

    cmd = ' '.join(['flirt', '-in', full_labels, '-ref', template,
                    '-applyxfm -init', xfm_matrix,
                    '-interp nearestneighbour -out', xfm_DKT31aseg])
    print(cmd); os.system(cmd)

