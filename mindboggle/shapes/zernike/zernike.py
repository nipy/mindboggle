#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(points, faces, n_moments=20):
    """
    Compute the Zernike moments of a surface patch of points and faces.
    
    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex
    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on a mesh
    n_moments : integer
        number of moments to compute

    Returns
    -------
    moments : numpy matrix
        Zernike moments

    Examples
    --------
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> vtk_file = '/drop/MB/data/arno/features/folds.vtk'
    >>> faces, u1,u2, points, u3,u4,u5,u6 = read_vtk(vtk_file)
    >>> #points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> #faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> n_moments = 20
    >>> zernike_moments(points, faces, n_moments)

    """
    import numpy as np
    from mindboggle.shapes.zernike.multiproc import MultiprocPipeline

    # Arguments should be numpy arrays:
    if isinstance(points, list):
        points = np.array(points)
    if isinstance(faces, list):
        faces = np.array(faces)

    pl = MultiprocPipeline()
    n_points = np.size(points)
    n_faces = np.size(faces)

    G = pl.geometric_moments_orig(points, faces, n_moments, n_faces, n_points)
    Z = pl.zernike(G, n_points)
    moments = pl.feature_extraction(Z, n_points)

    return moments


def zernike_moments_per_label(vtk_file, n_moments, exclude_labels=[-1]):
    """
    Compute the Zernike moments per labeled region in a file.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    n_moments : integer
        number of moments to compute
    exclude_labels : list of integers
        labels to be excluded

    Returns
    -------
    moments : numpy matrix
        Zernike moments
    label_list : list of integers
        list of unique labels for which moments are obtained

    Examples
    --------
    >>> import os
    >>> from mindboggle.shapes.laplace_beltrami import laplacian_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> n_moments = 20
    >>> zernike_moments_per_label(vtk_file, n_moments, exclude_labels=[-1]

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import remove_faces
    from mindboggle.shapes.zernike.zernike import zernike_moments

    # Read VTK surface mesh file:
    faces, u1,u2, points, u4, labels, u5,u6 = read_vtk(vtk_file)

    # Loop through labeled regions:
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    moments_lists = []
    for label in ulabels:
      #if label==22:

        # Determine the indices per label:
        label_indices = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.format(len(label_indices), label))

        # Remove background faces:
        select_faces = remove_faces(faces, label_indices)

        # Compute Zernike moments for the label:
        moments = zernike_moments(points, select_faces, n_moments)

        # Append to a list of lists of spectra:
        moments_lists.append(moments)
        label_list.append(label)

    return moments_lists, label_list
