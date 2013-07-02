#!/usr/bin/env python
"""
Rotate an image volume in every 90 degree orientation and view.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def rotate_view_repeat(input_file, ref_file):
    import os
    import numpy as np
    import nibabel as nb

    interp = 'nearestneighbour'

    view_command = 'itksnap'  # 'fslview'
    flirt_command = 'flirt'  # 'flirt'

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


    # Rotate 90 degrees in each direction
    print('Rotate 90 degrees in each direction')
    r = 'n'
    for irot, rot in enumerate(rotations):
        if r == 'y':
            break

        # Construct and save rotation matrix with translation about center of mass
        img = nb.load(input_file)
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
        matrix_file = 'rot' + str(irot) + '_' + os.path.basename(input_file) + '.txt'
        print('Save rotation matrix to file: ' + matrix_file)
        np.savetxt(matrix_file, matrix)

        # Apply rotation matrix and save and view output
        out_file = 'rot' + str(irot) + '_' + os.path.basename(input_file)
        print('Apply rotation matrix, save to file ' + out_file + ', and view')
        cmd = ' '.join([flirt_command, '-in', input_file,
                        '-ref', ref_file,
                        '-applyxfm -init', matrix_file,
                        '-interp', interp,
                        '-out', out_file])
        print(cmd); os.system(cmd)
        cmd = ' '.join([view_command, out_file, '&'])
        print(cmd); os.system(cmd)
