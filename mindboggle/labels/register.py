#!/usr/bin/env python
"""
Register VTK file points to template.


Authors:
- Arno Klein, 2013 (arno@mindboggle.info) http://binarybottle.com

Copyright 2013, Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def read_itk_transform(transform_file):
    """
    Read ITK transform file and output transform array.

    ..ITK affine transform file format ::

        #Insight Transform File V1.0
        #Transform 0
        Transform: MatrixOffsetTransformBase_double_3_3
        Parameters: 0.90768 0.043529 0.0128917 -0.0454455 0.868937 0.406098 \
                    0.0179439 -0.430013 0.783074 -0.794889 -18.3346 -3.14767
        FixedParameters: -0.60936 21.1593 10.6148

    Parameters
    ----------
    transform file : string
        name of ITK affine transform file

    Returns
    -------
    transform : numpy array
        4x4 affine transform matrix

    Examples
    --------
    >>> import os
    >>> from mindboggle.labels.register import read_itk_transform
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> transform_file = os.path.join(path, 'arno', 'mri',
    >>>                               't1weighted_brain.MNI152Affine.txt')
    >>> read_itk_transform(transform_file)
      array([[  9.07680000e-01,   4.35290000e-02,   1.28917000e-02, -6.09360000e-01],
             [ -4.54455000e-02,   8.68937000e-01,   4.06098000e-01, 2.11593000e+01],
             [  1.79439000e-02,  -4.30013000e-01,   7.83074000e-01, 1.06148000e+01],
             [ -7.94889000e-01,  -1.83346000e+01,  -3.14767000e+00, 1.00000000e+00]])

    """
    import numpy as np

    transform = np.eye(4)

    # Read ITK transform file
    fid = open(transform_file, 'r')
    affine_lines = fid.readlines()
    # Linear transform
    linear_transform = affine_lines[3]
    linear_transform = linear_transform.split()
    linear_transform = [np.float(x) for x in linear_transform[1::]]
    # Translation
    translation = affine_lines[4]
    translation = translation.split()
    translation = [np.float(x) for x in translation[1::]]

    # Write transform array
    transform[0:4,0:3] = np.reshape(linear_transform, (4,3))
    transform[0:3,3] = translation

    return transform

def transform_points(hemi, subject, transform,
                     subjects_path, atlas, atlas_string):
    """
    Transform the labels from a surface atlas via a template
    using FreeSurfer's mri_surf2surf (wrapped in NiPype).

    nipype.workflows.smri.freesurfer.utils.fs.SurfaceTransform
    wraps command ``mri_surf2surf``:

        "Transform a surface file from one subject to another via a spherical registration.
        Both the source and target subject must reside in your Subjects Directory,
        and they must have been processed with recon-all, unless you are transforming
        to one of the icosahedron meshes."

    Parameters
    ----------
    hemi : string
        hemisphere ['lh' or 'rh']
    subject : string
        subject corresponding to FreeSurfer subject directory
    transform : string
        name of FreeSurfer spherical surface registration transform file
    subjects_path : string
        name of FreeSurfer subjects directory
    atlas : string
        name of atlas
    atlas_string : string
        name of atlas labeling protocol

    """
#    from mindboggle.utils.io_vtk import read_points
#    points = read_points(file)
    from os import path, getcwd
    from nipype.interfaces.freesurfer import SurfaceTransform

    sxfm = SurfaceTransform()
    sxfm.inputs.hemi = hemi
    sxfm.inputs.target_subject = subject
    sxfm.inputs.source_subject = atlas

    # Source file
    sxfm.inputs.source_annot_file = path.join(subjects_path,
                                    atlas, 'label',
                                    hemi + '.' + atlas_string + '.annot')
    # Output annotation file
    output_file = path.join(getcwd(), hemi + '.' + atlas + '.' + atlas_string + \
                                      '_to_' + subject + '.annot')
    sxfm.inputs.out_file = output_file

    # Arguments: strings within registered files
    args = ['--srcsurfreg', transform,
            '--trgsurfreg', transform]
    sxfm.inputs.args = ' '.join(args)

    sxfm.run()

    return output_file
