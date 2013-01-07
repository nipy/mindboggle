#!/usr/bin/python
"""
Operations on volume matrices.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Flip axes
#------------------------------------------------------------------------------
def flip_axes(input_file, flipx=True, flipy=True, flipz=False,
              use_matrix=False, use_header=True):
    """
    Flip image volume data about x, y, and/or z axes
    without regard to the file header. The default is to flip
    x and y to rotate the image volume by 180 degrees.

    Parameters
    ----------
    input_file : image volume
    flipx : Boolean
        flip about x-axis?
    flipy : Boolean
        flip about y-axis?
    flipz : Boolean
        flip about z-axis?
    use_matrix : Boolean
        use input file's affine matrix?
    use_header : Boolean
        use input file's header?

    Returns
    -------
    out_file : string
        output image volume

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load image volume
    img = nb.load(input_file)
    dat = img.get_data()
    if use_matrix:
        mat = img.get_affine()
    if use_header:
        hdr = img.get_header()
    lenx, leny, lenz = np.shape(dat)
    dat_new = np.zeros((lenx, leny, lenz))

    # Flip x
    if flipx:
        for x in range(lenx):
            dat_new[lenx-1-x,:,:] = dat[x,:,:]

    # Flip y
    if flipy:
        for y in range(leny):
            dat_new[:,leny-1-y,:] = dat[:,y,:]

    # Flip z
    if flipz:
        for z in range(lenz):
            dat_new[:,:,lenz-1-z] = dat[:,:,z]

    # Save output
    out_file = 'reorient_' + os.path.basename(input_file)
    if use_matrix:
        if use_header:
            img = nb.Nifti1Image(dat_new, mat, hdr)
        else:
            img = nb.Nifti1Image(dat_new, mat)
    elif use_header:
        img = nb.Nifti1Image(dat_new, np.eye(4,4), hdr)
    else:
        img = nb.Nifti1Image(dat_new, np.eye(4,4))

    img.to_filename(out_file)

    return out_file

#------------------------------------------------------------------------------
# Flip axes
#------------------------------------------------------------------------------
def rotate90(input_file, rotations=[1,0,0],
             use_matrix=False, use_header=True):
    """
    Rotate image volume data 90 degrees about x, y, and/or z axes
    without regard to the file header.

    Parameters
    ----------
    input_file : image volume
    rotations : list of integers in {-1,0,1}
        rotate -90, 0, or +90 degrees about x-,y-,z-axis
    use_matrix : Boolean
        use input file's affine matrix?
    use_header : Boolean
        use input file's header?

    Returns
    -------
    out_file : string
        output image volume

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load image volume
    img = nb.load(input_file)
    dat = img.get_data()
    if use_matrix:
        mat = img.get_affine()
    if use_header:
        hdr = img.get_header()
    lenx, leny, lenz = np.shape(dat)
    dat_new = np.zeros((lenx, leny, lenz))

    minxyz = min([lenx, leny, lenz])

    # Rotate about x-axis
    for y in range(min([leny, lenz])):
        if rotations[0]==1:
            dat_new[0:minxyz, 0:minxyz, y] = dat[0:minxyz, y, 0:minxyz]
        elif rotations[0]==-1:
            dat_new[0:minxyz, 0:minxyz, leny-1-y] = dat[0:minxyz, y, 0:minxyz]

    # Rotate about y-axis
    for z in range(min([lenx, lenz])):
        if rotations[1]==1:
            dat_new[z, 0:minxyz, 0:minxyz] = dat[0:minxyz, 0:minxyz, z]
        elif rotations[1]==-1:
            dat_new[lenz-1-z, 0:minxyz, 0:minxyz] = dat[0:minxyz, 0:minxyz, z]

    # Rotate about z-axis
    for x in range(min([lenx, leny])):
        if rotations[1]==1:
            dat_new[0:minxyz, x, 0:minxyz] = dat[x, 0:minxyz, 0:minxyz]
        elif rotations[1]==-1:
            dat_new[0:minxyz, lenx-1-x, 0:minxyz] = dat[x, 0:minxyz, 0:minxyz]

    # Save output
    out_file = 'reorient_' + os.path.basename(input_file)
    if use_matrix:
        if use_header:
            img = nb.Nifti1Image(dat_new, mat, hdr)
        else:
            img = nb.Nifti1Image(dat_new, mat)
    elif use_header:
        img = nb.Nifti1Image(dat_new, np.eye(4,4), hdr)
    else:
        img = nb.Nifti1Image(dat_new, np.eye(4,4))

    img.to_filename(out_file)

    return out_file

#------------------------------------------------------------------------------
# Example
#------------------------------------------------------------------------------
if __name__ == "__main__":

    from mindboggle.utils.matrix_operations import rotate90

    input_file = 'PhStr_CU_20121130.nii.gz'
    rotate90(input_file, rotations=[1,0,0],
             use_matrix=False, use_header=False)


