#!/usr/bin/env python
"""
Utility functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def execute(cmd, type='os'):
    """
    Execute command by either subprocess.call or os.system.

    Parameters
    ----------
    cmd : sequence (string also permitted if type=='os')
        command with arguments
    type : string
        how to execute {os, subprocess}

    Examples
    --------
    >>> from mindboggle.utils.utils import execute
    >>> cmd = ['ls', '-l', '-a', '.']
    >>> type = 'subprocess'
    >>> execute(cmd, type)
    >>> type = 'os'
    >>> execute(cmd, type)
    >>> cmd = 'ls -l -a .'
    >>> execute(cmd)

    """
    from subprocess import call
    import sys

    if isinstance(cmd, str):
        print(cmd)
    else:
        print(' '.join(cmd))

    # Use subprocess.call:
    if type == 'subprocess':
        try:
            retcode = call(cmd)
            if retcode < 0:
                print >>sys.stderr, "Child terminated by signal", -retcode
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    # Use os.system:
    elif type == 'os':
        from os import system

        if isinstance(cmd, str):
            pass
        else:
            cmd = ' '.join(cmd)
        try:
            system(cmd)
        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    else:
        sys.exit('Select either "subprocess" or "os" for execution type.')


def list_strings(string1='', string2='', string3='', string4=''):
    """
    Put strings in a list.

    Parameters
    ----------
    string1 : string
    string2 : string
    string3 : string
    string4 : string

    Returns
    -------
    string_list : list of strings

    """

    string_list = []
    if string1:
        string_list.append(string1)
    if string2:
        string_list.append(string2)
    if string3:
        string_list.append(string3)
    if string4:
        string_list.append(string4)

    return string_list


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
