#!/usr/bin/env python
"""
Rotate an image volume in every 90 degree orientation and view.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

import os
import numpy as np
import nibabel as nb

# Inputs
dir = '/homedir/Data/EMBARC/Data'
#images_dir = os.path.join(dir, 'phantoms_dmri')
#file_list = os.path.join(dir, 'phantoms_dmri_1stvol.txt')
#file_list2 = os.path.join(dir, 'phantoms_dmri_FA.txt')
#ref_file = os.path.join(dir, 'phantoms_dmri_manufacturer/PhDifMFR_1stvol_CU_20120530.nii.gz')
#new_file_list = os.path.join(dir, 'phantoms_dmri_1stvol_rotated.txt')
#new_file_list2 = os.path.join(dir, 'phantoms_dmri_FA_rotated.txt')
images_dir = os.path.join(dir, 'phantoms_adni')
file_list = os.path.join(dir, 'phantoms_adni.txt')
file_list2 = os.path.join(dir, 'phantoms_adni.txt')
ref_file = os.path.join(dir, 'phantoms_adni/PhStr_CU_20120711.nii.gz')
new_file_list = os.path.join(dir, 'phantoms_adni_rotated.txt')
new_file_list2 = os.path.join(dir, 'phantoms_adni_rotated.txt')
interp = 'nearestneighbour'

view_command = 'itksnap'  # 'fslview'
flirt_command = 'flirt.fsl'  # 'flirt'

# Rotation matrices
rotations = np.array([[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
                    [[1,0,0,0],[0,0,-1,0],[0,1,0,0],[0,0,0,1]],
                    [[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]],
                    [[1,0,0,0],[0,0,1,0],[0,-1,0,0],[0,0,0,1]],
                    [[0,-1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]],
                    [[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]],
                    [[0,1,0,0],[1,0,0,0],[0,0,-1,0],[0,0,0,1]],
                    [[0,0,-1,0],[1,0,0,0],[0,-1,0,0],[0,0,0,1]],
                    [[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]],
                    [[-1,0,0,0],[0,0,-1,0],[0,-1,0,0],[0,0,0,1]],
                    [[-1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]],
                    [[-1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],
                    [[0,1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]],
                    [[0,0,1,0],[-1,0,0,0],[0,-1,0,0],[0,0,0,1]],
                    [[0,-1,0,0],[-1,0,0,0],[0,0,-1,0],[0,0,0,1]],
                    [[0,0,-1,0],[-1,0,0,0],[0,1,0,0],[0,0,0,1]],
                    [[0,0,-1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]],
                    [[0,1,0,0],[0,0,1,0],[1,0,0,0],[0,0,0,1]],
                    [[0,0,1,0],[0,-1,0,0],[1,0,0,0],[0,0,0,1]],
                    [[0,-1,0,0],[0,0,-1,0],[1,0,0,0],[0,0,0,1]],
                    [[0,0,-1,0],[0,-1,0,0],[-1,0,0,0],[0,0,0,1]],
                    [[0,-1,0,0],[0,0,1,0],[-1,0,0,0],[0,0,0,1]],
                    [[0,0,1,0],[0,1,0,0],[-1,0,0,0],[0,0,0,1]],
                    [[0,1,0,0],[0,0,-1,0],[-1,0,0,0],[0,0,0,1]]])

# Load file lists
fid_read = open(file_list, 'r')
lines = fid_read.readlines()
if len(lines) < 2:
    lines = lines[0].split()
fid_read2 = open(file_list2, 'r')
lines2 = fid_read2.readlines()
if len(lines2) < 2:
    lines2 = lines2[0].split()

# Save new file list
fid_write = open(new_file_list, 'w')

# Process each file
for iline, line in enumerate(lines):
#  if 'CU_20121130' in line:
    file = line.strip()
    full_path = os.path.join(images_dir, file)

    # Rotate 90 degrees in each direction
    print('Rotate 90 degrees in each direction')
    r = 'n'
    for irot, rot in enumerate(rotations):
        if r == 'y':
            break

        # Construct and save rotation matrix with translation about center of mass
        img = nb.load(full_path)
        hdr = img.get_header()
        pixdims = hdr.get_zooms()
        dat = img.get_data()
        i,j,k = np.where(dat > 0.1)
        center = np.array((sum(i), sum(j), sum(k))) / len(i)
        shift = pixdims * center
        trans1 = np.eye(4,4)
        trans2 = np.eye(4,4)
        trans1[0:3,-1] = -shift
        trans2[0:3,-1] = shift
        matrix = np.dot(trans2, np.dot(rot, trans1))
        matrix_file = 'rot' + str(irot) + '_' + os.path.basename(file) + '.txt'
        print('Save rotation matrix to file: ' + matrix_file)
        np.savetxt(matrix_file, matrix)

        # Apply rotation matrix and save and view output
        out_file = 'rot' + str(irot) + '_' + os.path.basename(file)
        print('Apply rotation matrix, save to file ' + out_file + ', and view')
        cmd = ' '.join([flirt_command, '-in', os.path.join(images_dir, file),
                        '-ref', ref_file,
                        '-applyxfm -init', matrix_file,
                        '-interp', interp,
                        '-out', out_file])
        print(cmd); os.system(cmd)
        cmd = ' '.join([view_command, out_file, '&'])
        print(cmd); os.system(cmd)

        # Skip to the next file if this is the correct rotation
        r = raw_input('Did this result in the correct orientation? ("y" or press Return)\n\n\t')
        if r == 'y':
            file2 = lines2[iline].strip()
            out_file2 = 'rot' + str(irot) + '_' + os.path.basename(file2)
            print('Apply rotation matrix, save to file ' + out_file2)
            cmd = ' '.join([flirt_command, '-in', os.path.join(images_dir, file2),
                            '-ref', ref_file,
                            '-applyxfm -init', matrix_file,
                            '-interp', interp,
                            '-out', out_file2])
            print(cmd); os.system(cmd)

            print('The correctly initialized file is ' + out_file2)
            fid_write = open(new_file_list, 'a')
            fid_write.write(out_file + '\n')
            fid_write.close()
            fid_write2 = open(new_file_list2, 'a')
            fid_write2.write(out_file2 + '\n')
            fid_write2.close()
            break
