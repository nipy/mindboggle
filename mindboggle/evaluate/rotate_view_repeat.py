
import os
import numpy as np
import nibabel as nb
from mindboggle.utils.matrix_operations import make_voxels_isometric, rotate90

# Inputs
dir = '/drop/EMBARC/Data/Phantom_DTI_Updated20121214/'
file_list = os.path.join(dir, 'Phantom_DTI_1stVol.txt')
file_list2 = os.path.join(dir, 'Phantom_DTI_FA.txt')
ref_file = '/drop/EMBARC/Data/Phantom_DTI_Manufacturer_Updated20130108/PhDifMFR_1stvol_CU_20120530.nii.gz'
new_file_list = os.path.join(dir, 'Phantom_DTI_1stVol_rotated.txt')
new_file_list2 = os.path.join(dir, 'Phantom_DTI_FA_rotated.txt')
interp = 'nearestneighbour'

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
#for iline, line in enumerate(lines[0::]):
#    iline=0
    file = line.strip()
    full_path = os.path.join(dir, file)
    #file = os.path.basename(line)
    #dir = '/'.join(line.split('/')[0:-1])

    # Rotate 90 degrees in each direction
    print('Rotate 90 degrees in each direction')
    r = 'n'
    for irot, rot in enumerate(rotations):
        if r == 'y':
            break

        # Construct and save rotation matrix with translation about center of mass
        img = nb.load(full_path)
        #dims = img.get_shape()
        #halfdims = np.array(dims) / 2
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
        cmd = ' '.join(['flirt -in', os.path.join(dir, file),
                        '-ref', ref_file,
                        '-applyxfm -init', matrix_file,
                        '-interp', interp,
                        '-out', out_file])
        print(cmd); os.system(cmd)
        cmd = 'fslview ' + out_file + ' &'
        print(cmd); os.system(cmd)

        # Skip to the next file if this is the correct rotation
        r = raw_input('Did this result in the correct orientation? ("y" or press Return)\n\n\t')
        if r == 'y':
            file2 = lines2[iline].strip()
            out_file2 = 'rot' + str(irot) + '_' + os.path.basename(file2)
            print('Apply rotation matrix, save to file ' + out_file2)
            cmd = ' '.join(['flirt -in', os.path.join(dir, file2),
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
