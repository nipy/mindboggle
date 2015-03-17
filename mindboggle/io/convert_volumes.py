#!/usr/bin/env python
"""
Functions for reading and writing nifti volume files.


Authors:
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def convert2nii(input_file, reference_file, output_file='', interp='continuous'):
    """
    Convert volume from the input file format to the output file format.

    If output_file is empty, reslice to nifti format using nibabel and
    scipy.ndimage.affine_transform, after nilearn.image.resample_img::

        from nilearn.image import resample_img
        resliced = resample_img(input_file, target_affine=xfm2,
                                target_shape=dim2,
                                interpolation=interp).get_data()

    Example use: Convert FreeSurfer 'unconformed' .mgz file to nifti.

    Parameters
    ----------
    input_file : string
        input file name
    reference_file : string
        target file name
    output_file : string
        name of output file
    interp : string
        interpolation method: 'continuous' (default) or 'nearest'

    Returns
    -------
    output_file : string
        name of output file

    Examples
    --------
    >>> import os
    >>> from mindboggle.io.nifti_io import convert2nii
    >>> subject = 'Twins-2-1'
    >>> input_file = os.path.join(os.environ['SUBJECTS_DIR'],subject,'mri','aparc+aseg.mgz')
    >>> reference_file = os.path.join(os.environ['SUBJECTS_DIR'],subject,'mri','orig','001.mgz')
    >>> output_file = ''
    >>> interp = 'nearest'
    >>> output_file = convert2nii(input_file, reference_file, output_file, interp)
    >>> #from mindboggle.io.plots import plot_volumes
    >>> #command = '/Applications/ITK-SNAP.app/Contents/MacOS/InsightSNAP'
    >>> #plot_volumes('orig.mgz.nii.gz', command=command)

    """
    import os
    import numpy as np
    import nibabel as nb
    from scipy import ndimage, linalg

    if not os.path.exists(input_file):
        raise(IOError("Input file " + input_file + " not found"))
    if not os.path.exists(reference_file):
        raise(IOError("Reference file " + reference_file + " not found."))
    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(input_file) + '.nii.gz')
    #-------------------------------------------------------------------------
    # Load reference image:
    #-------------------------------------------------------------------------
    vol2 = nb.load(reference_file)
    xfm2 = vol2.get_affine()
    dim2 = vol2.shape
    #-------------------------------------------------------------------------
    # Resample the source image according to the reference image:
    #-------------------------------------------------------------------------
    vol1 = nb.load(input_file)
    dat1 = vol1.get_data()
    xfm1 = vol1.get_affine()
    if np.all(xfm2 == xfm1):
        transform_affine = np.eye(4)
    else:
        transform_affine = np.dot(linalg.inv(xfm1), xfm2)
    A = transform_affine[0:3, 0:3]
    b = transform_affine[0:3, 3]
    A_inv = linalg.inv(A)
    # If A is diagonal, affine_transform uses a better algorithm.
    if np.all(np.diag(np.diag(A)) == A):
        A = np.diag(A)
    else:
        b = np.dot(A, b)

    # order of the spline interpolation:
    if interp == 'nearest':
        interpolation_order = 0
    else:
        interpolation_order = 3
    resliced = ndimage.affine_transform(dat1, A,
                         offset=np.dot(A_inv, b),
                         output_shape=dim2,
                         order=interpolation_order)

    #-------------------------------------------------------------------------
    # Save the image with the reference affine transform:
    #-------------------------------------------------------------------------
    img = nb.Nifti1Image(resliced, xfm2)
    img.to_filename(output_file)

    return output_file


def xyz2nii(input_xyz_file, output_nii_file='', origin=[], pad=10):
    """
    Convert [x,y,z] coordinate file to nifti (nii.gz) volume file.

    Parameters
    ----------
    input_xyz_file : string
        input [x,y,z] coordinate text file
    output_nii_file : string
        output nifti (nii.gz) volume file
    origin : list of floats
        [x,y,z] coordinates for origin
    pad : integer
        number of voxels to pad input coordinates in x, y, and z directions

    Returns
    -------
    output_nii_file : string
        output nifti (nii.gz) volume file

    Examples
    --------
    >>> from mindboggle.io.nifti_io import xyz2nii
    >>> input_xyz_file = '/Users/arno/Dropbox/MSSM/Nebojsa/face.xyz.txt'
    >>> origin = []
    >>> pad = 10
    >>> output_nii_file = ''
    >>> xyz2nii(input_xyz_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load coordinates and scalars:
    XYZscalars = np.loadtxt(input_xyz_file)
    XYZ = np.round(XYZscalars[:, 0:3])
    #scalars = XYZscalars[:, 3::]

    if origin:
        XYZ -= origin

    XYZ += np.abs(np.min(XYZ, axis=0)) + [pad, pad, pad]
    XYZ = np.round(XYZ)
    dims = np.max(XYZ, axis=0) + [pad, pad, pad]
    data = np.zeros(dims)

    # Loop through rows or array and write 1s in image volume:
    for irow, xyz in enumerate(XYZ):
        data[xyz[0], xyz[1], xyz[2]] = 1

    # Write output image volume:
    if not output_nii_file:
        output_nii_file = os.path.join(os.getcwd(), 'xyz.nii.gz')
    img = nb.Nifti1Image(data, affine=np.eye(4,4))
    img.to_filename(output_nii_file)

    return output_nii_file
