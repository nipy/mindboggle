#!/usr/bin/python
"""
Count landmarks.

(c) 2011, @rno klein
"""

import os
from os.path import exists
from subprocess import call
import numpy as np
from nibabel import load

verbose = 1
dim = 3
files  = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']

out_path    = '/hd2/Archive/registration_evaluation_2011_output/'
#in_path    = os.path.join( out_path, 'Transformed_Landmarks/')
in_path     = '/hd2/Brains/CUMC12/'
results_dir = os.path.join(out_path, 'Results/')
ext         = '.nii.gz'

landmark_type = 'fundi_gang_li'
landmark_dir = in_path+'Landmarks/fundi_gang_li_binary/'
landmark_type = 'fundi_brain_visa'
landmark_dir = in_path+'Landmarks/fundi_brain_visa_binary/'
landmark_type = 'ribbons_brain_visa'
landmark_dir = in_path+'Landmarks/ribbons_brain_visa_binary/'
landmark_type = 'pits_forrest_bao'
landmark_dir = in_path+'Landmarks/pits_forrest_bao_binary/'
landmark_type = 'pits_yrjo_hame'
landmark_dir = in_path+'Landmarks/pits_yrjo_hame_binary/'
landmark_type = 'pits_kiho_im'
landmark_dir = in_path+'Landmarks/pits_kiho_im_binary/'
landmark_type = 'fundi_forrest_bao'
landmark_dir = in_path+'Landmarks/fundi_forrest_bao_binary/'

print(landmark_type)
results_file = results_dir+'count_'+landmark_type+'.txt'
count_voxels = []
count_pieces = []
temp_file = results_dir+'temp_count_'+landmark_type+'.txt'
f = open(results_file, 'w');
f.write('Count: voxels, components \n')
for file in files:
    #source_file = landmark_dir+file+'_to_template_'+landmark_type+ext
    source_file = landmark_dir+file+ext
    out_file = landmark_dir+'connected_'+file+ext
    if os.path.exists(source_file):

        I = load(source_file)
        I = I.get_data()
        count1 = int(sum(np.ravel(I)))

        args = " ".join(['c3d', source_file, '-connected-components -o', out_file])
        p = call(args, shell="True")
        args = " ".join(['c3d', out_file, out_file, '-label-statistics > ', temp_file])
        p = call(args, shell="True")
        ftemp = open(temp_file, 'r')
        count2 = int(len(ftemp.readlines()) - 2)
        
        print_out = " ".join([file+':', str(count1), str(count2)])
        print(print_out)
        f.close()
        f = open(results_file, 'a')
        f.write(print_out + '\n')
        count_voxels.append(count1)
        count_pieces.append(count2)
    else:
        if not os.path.exists(source_file):
            raise NameError('Check input file ' + source_file)

avg_count1 = str(np.mean(count_voxels))
std_count1 = str(np.std(count_voxels))
min_count1 = str(int(np.min(count_voxels)))
max_count1 = str(int(np.max(count_voxels)))
avg_count2 = str(np.mean(count_pieces))
std_count2 = str(np.std(count_pieces))
min_count2 = str(int(np.min(count_pieces)))
max_count2 = str(int(np.max(count_pieces)))

print_out1 = " ".join(['Number of voxels (mean, std, min, max):', \
                      avg_count1, std_count1, min_count1, max_count1])
print_out2 = " ".join(['Number of pieces (mean, std, min, max):', \
                      avg_count2, std_count2, min_count2, max_count2])
print(print_out1)
print(print_out2)
f.close()
f = open(results_file, 'a');
f.write(print_out1 + '\n' + print_out2)
f.close()
