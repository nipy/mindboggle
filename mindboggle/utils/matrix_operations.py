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
    without regard to the file header.  Assumes isometric voxels.

    Parameters
    ----------
    input_file : image volume
    rotations : list of integers in {-1,0,1}
        rotate -90, 0, and/or +90 degrees about x-, y-, or z-axis (in order)
    use_matrix : Boolean
        use input file's affine matrix?
    use_header : Boolean
        use input file's header?

    Returns
    -------
    out_file : string
        output image volume

    Examples
    --------
    >>> # Example without proper inputs specified:
    >>> # Make voxels isometric:
    >>> import os
    >>> d = os.listdir('.')
    >>> for f in d:
    >>>     if '1stvol' in f:
    >>>         c=' '.join(['c3d',f,'-resample 128x128x128','-o isometric_'+f])
    >>>         print(c); os.system(c)
    >>> # Rotate 90 degrees about x-axis:
    >>> from mindboggle.utils.matrix_operations import rotate90
    >>> rotate90(input_file, rotations=[1,0,0], use_matrix=False, use_header=True)

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
    if lenx != leny or lenx != lenz or leny != lenz:
        os.error('Input file should have isometric voxels.')

    dat_new = np.zeros((lenx, leny, lenz))

    # Rotate about x-axis
    if rotations[0] != 0:
        for y in range(lenx):
            if rotations[0]==1:
                dat_new[:, :, y] = dat[:, y, :]
            elif rotations[0]==-1:
                dat_new[:, :, leny-1-y] = dat[:, y, :]

    # Rotate about y-axis
    if rotations[1] != 0:
        for z in range(lenx):
            if rotations[1]==1:
                dat_new[z, :, :] = dat[:, :, z]
            elif rotations[1]==-1:
                dat_new[lenz-1-z, :, :] = dat[:, :, z]

    # Rotate about z-axis
    if rotations[2] != 0:
        for x in range(lenx):
            if rotations[1]==1:
                dat_new[:, x, :] = dat[x, :, :]
            elif rotations[1]==-1:
                dat_new[:, lenx-1-x, :] = dat[x, :, :]

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

    import os
    from mindboggle.utils.matrix_operations import rotate90

    input_dir = '/drop/EMBARC/Data/Phantom_DTI_Updated20121214'
    input_files = ['isometric_PhDif_1stvol_CU_20120607.nii.gz',
                   'isometric_PhDif_1stvol_CU_20120711.nii.gz',
                   'isometric_PhDif_1stvol_CU_20120915.nii.gz',
                   'isometric_PhDif_1stvol_CU_20121206.nii.gz',
                   'isometric_PhDif_1stvol_UM_20120803.nii.gz',
                   'isometric_PhDif_1stvol_UM_20121105.nii.gz',
                   'isometric_PhDif_1stvol_UM_20121204.nii.gz']

    for input_file in input_files:
        input_file = os.path.join(input_dir, input_file)
        rotate90(input_file, rotations=[1,0,0],
                 use_matrix=False, use_header=False)


