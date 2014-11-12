#!/usr/bin/env python
"""
Functions for reading and writing nifti volume files.


Authors:
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def convert2nii(input_file, reference_file, output_file='', interp=1):
    """
    Convert volume from the input file format to the output file format.
    If output_file is empty, reslice to nifti format using nibabel and
    scipy.ndimage.

    Example use: Convert FreeSurfer 'unconformed' .mgz file to nifti.

    Parameters
    ----------
    input_file : string
        input file name
    reference_file : string
        file in original space
    output_file : string
        name of output file
    interp : integer
        interpolation method: 0=nearest, 1=trilin, >1=spline

    Returns
    -------
    output_file : string
        name of output file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.freesurfer import convert_mgh_to_native_nifti
    >>> from mindboggle.utils.plots import plot_volumes
    >>> subject = 'bert'
    >>> input_file = os.path.join(os.environ['SUBJECTS_DIR'],subject,'mri','orig','001.mgz')
    >>> reference_file = input_file
    >>> convert_mgh_to_native_nifti(input_file, reference_file)
    >>> command = '/Applications/ITK-SNAP.app/Contents/MacOS/InsightSNAP'
    >>> plot_volumes('001.nii.gz', command=command)

    """
    import os
    import numpy as np
    import nibabel as nb
    from scipy.ndimage.interpolation import affine_transform

    print("Convert volume from FreeSurfer 'unconformed' to original space...")

    if not os.path.exists(input_file):
        raise(IOError("Input file " + input_file + " not found"))
    if not os.path.exists(reference_file):
        raise(IOError("Reference file " + reference_file + " not found."))
    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   os.path.basename(input_file) + 'nii.gz')

    # Load image volume
    vol = nb.load(input_file)
    dat = vol.get_data()
    dim = vol.shape
    xfm = vol.get_affine()
    xfm2 = np.eye(3,3)

    resliced = affine_transform(dat, xfm2, output_shape=dim, order=interp)

    # Save the image with the same affine transform:
    img = nb.Nifti1Image(resliced, xfm)
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
    >>> from mindboggle.utils.io_table import xyz2nii
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
