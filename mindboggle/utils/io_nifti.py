#!/usr/bin/env python
"""
Functions for reading and writing nifti image volumes.


Authors:
    - Arno Klein, 2013-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


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
