#!/usr/bin/python
"""
Run registration commands
(c) 2011, @rno klein
"""
import sys
import os
import os.path

#------------
# Directories
#------------
source_files       = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
target_files       = source_files
ANTSPATH           = '/software/ANTS-build/'
out_path           = '/hd2/data/Archive/registration_evaluation_2011_output/'
xfm_dir            = os.path.join( out_path, 'Transforms/')
xfm_atlas_dir      = os.path.join( out_path, 'Transformed_Atlases/')
atlas_dir          = '/hd2/data/Brains/CUMC12/Atlases/'
brain_dir          = '/hd2/data/Brains/CUMC12/Brains/'
atlas_append       = '.nii.gz'
brain_append       = '.nii.gz'
file_append        = '.nii.gz'
ants_template      = '/hd2/data/Brains/CUMC12/CUMC12_template.nii.gz'
xfm_via_template   = 1

#------------------------------
# Register subjects to template
#------------------------------
for file in source_files:

    cmd = ANTSPATH + 'ANTS 3 -m CC[' + ants_template + ',' + \
          brains_dir+file+brains_append+',1,3] -o ' + xfm_dir+'/'+file+file_append + \
          ' -r Gauss[3,0] -t SyN[0.5]  -i 30x100x10 --use-Histogram-Matching'
    print cmd; #os.system(cmd)

#----------------------
# Pairwise registration
#----------------------
"""
for file in source_files:
    for file2 in target_files:
    
        cmd = ANTSPATH + 'WarpImageMultiTransform 3 ' + atlas_dir+file+atlas_append + ' ' + \
              xfm_atlas_dir+reg_method+'/'+file+'_to_'+file2+file_append + \
              ' -R ' + atlas_dir+file2+atlas_append + \
              ' -i ' + xfm_dir+file2+'Affine.txt ' + \
              xfm_dir+file2+'InverseWarp.nii.gz ' + \
              xfm_dir+file+'Warp.nii.gz ' + \
              xfm_dir+file+'Affine.txt --use-NN '
        print cmd; os.system(cmd)
"""      
      
