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
    import os
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

def apply_affine_transform(transform_file, vtk_file):
    """
    Transform coordinates using an affine matrix.

    Parameters
    ----------
    transform file : string
        name of ITK affine transform file
    vtk_file : string
        name of VTK file containing point coordinate data

    Returns
    -------
    affined_points : list of lists of floats
        transformed coordinates
    output_file : string
        name of VTK file containing transformed point data

    Examples
    --------
    >>> import os
    >>> from mindboggle.labels.register import apply_affine_transform
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> transform_file = os.path.join(path, 'arno', 'mri',
    >>>                               't1weighted_brain.MNI152Affine.txt')
    >>> vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> apply_affine_transform(transform_file, vtk_file)
    >>> # View
    >>> plot_vtk('transformed_lh.pial.depth.vtk')


    """
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk, write_vtk
    from mindboggle.labels.register import read_itk_transform

    # Read ITK affine transform file
    transform = read_itk_transform(transform_file)

    # Read VTK file
    faces, lines, indices, points, npoints, scalars, name = read_vtk(vtk_file)

    # Transform points
    affined_points = transform * np.transpose(points)

    # Output transformed VTK file
    output_file = os.path.join(os.getcwd(), 'affined_' + os.path.basename(vtk_file))

    # Write VTK file
    write_vtk(output_file, affined_points, indices, lines, faces, scalars, name)

    return affined_points, output_file


if __name__ == "__main__" :
    import os
    from mindboggle.labels.register import apply_affine_transform
    from mindboggle.utils.mesh import plot_vtk
    path = os.environ['MINDBOGGLE_DATA']
    transform_file = os.path.join(path, 'arno', 'mri',
                                  't1weighted_brain.MNI152Affine.txt')
    vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    apply_affine_transform(transform_file, vtk_file)
    # View
    plot_vtk('transformed_lh.pial.depth.vtk')
