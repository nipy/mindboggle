#!/usr/bin/python
"""
Copy FreeSurfer and manually labeled files to populate
Mindboggle's data directory or the Mindboggle-101 distribution.

Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""
import os

run_this = 5 # 1 = Mindboggle-101-volumes
             # 2 = Mindboggle-101-surfaces
             # 3 = Mindboggle-101-FS: FreeSurfer files / Mindboggle + evaluation
             # 4 = atlas registration in Mindboggle
             # 5 = Mindboggle-101-volumes in MNI152 space

#------------------------------------------------
# Subjects file and source and target directories
#------------------------------------------------
list_file = '/projects/Mindboggle/mindboggle/mindboggle/info/atlases101.txt'
srcs_path = '/hd2/Lab/Brains/Mindboggle101/subjects'
maindir = '/projects/Mindboggle/Mindboggle-101'
if run_this == 1:
    tgts_path = os.path.join(maindir, 'volumes')
elif run_this == 2:
    tgts_path = os.path.join(maindir, 'surfaces')
elif run_this == 3:
    tgts_path = os.path.join(maindir, 'FS')
elif run_this == 4:
    tgts_path = os.path.join(maindir, 'mindboggle_data')
elif run_this == 5:
    tgts_path = os.path.join(maindir, 'volumes_in_MNI152')
os.system('mkdir ' + maindir)
os.system(tgts_path)

#--------------
# Files to copy
#--------------
# Mindboggle-101-volumes distribution
if run_this == 1:
    src_dirs = ['mri']
    tgt_dirs = []
    copy_files = [['t1weighted.nii.gz','t1weighted_brain.nii.gz',
                   'labels.DKT25.manual.nii.gz','labels.DKT31.manual.nii.gz']]
# Mindboggle-101-surfaces distribution
if run_this == 2:
    src_dirs = ['label']
    tgt_dirs = []
    copy_files = [['lh.labels.DKT25.manual.vtk','rh.labels.DKT25.manual.vtk',
                   'lh.labels.DKT31.manual.vtk','rh.labels.DKT31.manual.vtk']]
# Mindboggle-101 distribution FreeSurfer supplement;
# everything needed to run Mindboggle as atlas or as target, and evaluate:
elif run_this == 3:
    src_dirs = ['mri','surf','label']
    tgt_dirs = src_dirs
    copy_files = [['aseg.mgz','ribbon.mgz'],
                  ['lh.pial','lh.sphere','lh.inflated','lh.white',
                   'rh.pial','rh.sphere','rh.inflated','rh.white'],
                  ['lh.labels.DKT25.manual.annot','rh.labels.DKT25.manual.annot']]
# Everything needed to run just Mindboggle atlas registration (no evaluation):
elif run_this == 4:
    src_dirs = ['surf','label']
    tgt_dirs = src_dirs
    copy_files = [['lh.sphere','rh.sphere'],
                  ['lh.labels.DKT25.manual.annot','rh.labels.DKT25.manual.annot']]
# Mindboggle-101-volumes in MNI152 space
elif run_this == 5:
    src_dirs = ['mri']
    tgt_dirs = []
    copy_files = [['t1weighted_brain.MNI152.nii.gz','t1weighted_brain.MNI152.mat',
                   'labels.DKT25.manual.MNI152.nii.gz','labels.DKT31.manual.MNI152.nii.gz']]

#---------------------------------------------
# Copy files from source to target directories
#---------------------------------------------
if not os.path.exists(tgts_path):
    os.system('mkdir ' + tgts_path)

# For each subject in the subjects file
fid = open(list_file, 'r')
subjects = fid.readlines()
subjects = [''.join(x.split()) for x in subjects]
#subjects = ['Afterthought-1']
#subjects = ['OASIS-TRT-20-11']
for subject in subjects:
    print(subject)
    tgt_path = os.path.join(tgts_path, subject)
    if not os.path.exists(tgt_path):
        os.system('mkdir ' + tgt_path)

    # For each directory to copy from
    for i, src_dir in enumerate(src_dirs):
        src_path = os.path.join(srcs_path, subject, src_dir)
        if len(tgt_dirs):
            tgt_path = os.path.join(tgts_path, subject, tgt_dirs[i])
            if not os.path.exists(tgt_path):
                os.system('mkdir ' + tgt_path)

        # For each file to copy
        for copy_file in copy_files[i]:
            os.system(' '.join(['cp', os.path.join(src_path, copy_file), tgt_path]))
        #print(' '.join(['rm', os.path.join(src_path, 't1weighted.MNI152Affine.txt')]))
        #os.system(' '.join(['rm', os.path.join(src_path, 't1weighted.MNI152.nii.gz')]))

