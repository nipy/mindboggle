#!/usr/bin/env python
"""
Operations on volume matrices.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def make_voxels_isometric(input_file):
    """
    Pad smaller dimensions to reach maximum dimension (max_dim).
    Resample so that the x, y, and z dimensions are equal to the max_dim.

    Parameters
    ----------
    input_file : string
        nibabel-readable image volume file

    Returns
    -------
    out_file : string
        output isometric image volume

    Examples
    --------
    >>> # Example without proper inputs specified:
    >>> import os
    >>> from mindboggle.utils.matrix import make_voxels_isometric
    >>> d = os.listdir('.')
    >>> for f in d:
    >>>     make_voxels_isometric(f)

    """
    import os
    import nibabel as nb
    from mindboggle.utils.utils import execute

    # Load image volume
    img = nb.load(input_file)
    lenx, leny, lenz = img.shape
    if lenx != leny or lenx != lenz or leny != lenz:

        # Pad smaller dimensions
        max_dim = max([lenx, leny, lenz])
        padx = max_dim - lenx
        pady = max_dim - leny
        padz = max_dim - lenz

        # Save padded output
        out_file = 'isometric_' + os.path.basename(input_file)
        pad_dims = '{0}x{1}x{2}vox'.format(padx, pady, padz)
        cmd = ['c3d', input_file, '-pad 0x0x0vox', pad_dims, '-o ', out_file]
        execute(cmd)

        # Resample output
        max_dims = '{0}x{1}x{2}'.format(max_dim, max_dim, max_dim)
        cmd = ['c3d', out_file, '-resample', max_dims, '-o ', out_file]
        execute(cmd)

    else:
        out_file = input_file

    return out_file


def flip_axes(input_file, flipx=True, flipy=True, flipz=False,
              use_matrix=False, use_header=True):
    """
    Flip image volume data about x, y, and/or z axes
    without regard to the file header. The default is to flip
    x and y to rotate the image volume by 180 degrees.

    Parameters
    ----------
    input_file : string
        nibabel-readable image volume file
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


def rotate90(input_file, rotations=[1,0,0],
             use_matrix=True, use_header=True):
    """
    Pad and rotate image volume data 90 degrees about the x, y, or z axis
    without regard to the file header.

    Note ::
        Assumes isometric voxels and equal dimensions along each axis.

    Parameters
    ----------
    input_file : string
        nibabel-readable image volume file
    rotations : list of integers in {0,1}
        rotate 0, or +90 degrees about x-, y-, or z-axis
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
    >>> import os
    >>> from mindboggle.utils.matrix import rotate90
    >>> d = os.listdir('.')
    >>> for f in d:
    >>>     # Rotate 90 degrees about x-axis:
    >>>     rotate90(f, rotations=[1,0,0], use_matrix=False, use_header=True)

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.utils.utils import execute

    # Load image volume
    img = nb.load(input_file)
    dat = img.get_data()
    if use_matrix:
        mat = img.get_affine()
    if use_header:
        hdr = img.get_header()
        pixdims = hdr.get_zooms()
        if pixdims[0] != pixdims[1] or \
           pixdims[0] != pixdims[2] or \
           pixdims[1] != pixdims[2]:
            os.error('This function assumes isometric voxels.')
    else:
        pixdims = [1,1,1]
    dims = img.get_shape()
    if dims[0] != dims[1] or \
       dims[0] != dims[2] or \
       dims[1] != dims[2]:
        os.error('This function assumes same dimensions along each axis.')

    dat_new = np.zeros(dims)
    lenx = dims[0]

    # Rotate about x-axis
    if rotations[0] != 0:
        print('Rotate about x-axis')
        pixdims = '{0} {1} {2}'.format(pixdims[0], pixdims[2], pixdims[1])
        for y in range(lenx):
            if rotations[0]==1:
                dat_new[:, :, y] = dat[:, y, :]
            #elif rotations[0]==-1:
            #    dat_new[:, :, lenx-1-y] = dat[:, y, :]

    # Rotate about y-axis
    elif rotations[1] != 0:
        print('Rotate about y-axis')
        pixdims = '{0} {1} {2}'.format(pixdims[2], pixdims[1], pixdims[0])
        for z in range(lenx):
            if rotations[1]==1:
                dat_new[z, :, :] = dat[:, :, z]
            #elif rotations[1]==-1:
            #    dat_new[lenx-1-z, :, :] = dat[:, :, z]

    # Rotate about z-axis
    elif rotations[2] != 0:
        print('Rotate about z-axis')
        pixdims = '{0} {1} {2}'.format(pixdims[1], pixdims[0], pixdims[2])
        for x in range(lenx):
            if rotations[2]==1:
                dat_new[:, x, :] = dat[x, :, :]
            #elif rotations[2]==-1:
            #    dat_new[:, lenx-1-x, :] = dat[x, :, :]

    # Save output
    out_file = 'rot'+''.join([str(x) for x in rotations])\
                 + '_' + os.path.basename(input_file)
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

    # Swap pixel dimensions
    print('Swap pixel dimensions')
    cmd = ['fslchpixdim', out_file, pixdims]
    execute(cmd)

    return out_file


def crop_to_match_volume(source, target, output=''):
    """
    Crop source volume to match target volume.

    Center source image in target dimensions; save with target header.

    Parameters
    ----------
    source : string
        nibabel-readable (e.g., nifti) file
    target : string
        nibabel-readable (e.g., nifti) file
    output : string
        nibabel-readable (e.g., nifti) file

    Returns
    -------
    output : string
        nibabel-readable (e.g., nifti) file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.matrix import crop_to_match_volume
    >>> from mindboggle.utils.plots import plot_volumes
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> source = os.path.join(data_path, 'arno', 'labels', 'labels.DKT31.manual.nii.gz')
    >>> target = os.path.join(data_path, 'atlases', 'MNI152_T1_1mm_brain.nii.gz')
    >>> output = ''
    >>> ignore_labels = [0]
    >>> output = crop_source_to_match_target_volume(source, target, output)
    >>> # View
    >>> plot_volumes(output)

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load source and target image volumes:
    vol_source = nb.load(source)
    vol_target = nb.load(target)
    source_shape = vol_source.shape
    target_shape = vol_target.shape
    xfm = vol_target.get_affine()
    data_source = vol_source.get_data()

    # Crop source to target dimensions:
    shifts = [np.floor((source_shape[0] - target_shape[0]) / 2),
              np.floor((source_shape[1] - target_shape[1]) / 2),
              np.floor((source_shape[2] - target_shape[2]) / 2)]
    data_source_crop = data_source[shifts[0]:target_shape[0]+shifts[0],
                                   shifts[1]:target_shape[1]+shifts[1],
                                   shifts[2]:target_shape[2]+shifts[2]]

    # Save cropped source file:
    if not output:
        output = os.path.join(os.getcwd(),
                              os.path.basename(source).split('.')[0] +
                              '_to_' +
                              os.path.basename(target).split('.')[0] +
                              '.nii.gz')
    img = nb.Nifti1Image(data_source_crop, xfm)
    img.to_filename(output)

    return output


#-----------------------------------------------------------------------------
# Example
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    from mindboggle.utils.matrix import rotate90

    input_dir = '/drop/phantoms_dmri'
    input_files = ['isometric_PhDif_1stvol_CU_20120607.nii.gz',
                   'isometric_PhDif_1stvol_UM_20121204.nii.gz']
    input_dir = '/drop/phantoms_adni'
    input_files = ['PhStr_CU_20121130.nii.gz',
                   'PhStr_raw_CU_20130204.nii.gz']

    for input_file in input_files:
        input_file = os.path.join(input_dir, input_file)
        rotate90(input_file, rotations=[1,0,0],
                 use_matrix=False, use_header=False)


