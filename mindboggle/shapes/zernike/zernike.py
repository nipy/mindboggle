#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(points, faces, order=20):
    """
    Compute the Zernike moments of a surface patch of points and faces.
    
    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex
    faces : list of lists of 3 integers
        each list contains indices to vertices that form a triangle on a mesh
    order : integer
        order of the moments being calculated

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors

    Examples
    --------
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> order = 3
    >>> zernike_moments(points, faces, order)
    [0.07957747154594767,
     0.045583723371562586,
     0.12969541097557186,
     0.09450348267576704,
     0.1377604175344217,
     0.12763442544037504]
    >>> # Moments for label 22 (postcentral) in Twins-2-1
    >>> # (after running explode_scalars() with reindex=True):
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> vtk_file = '/drop/MB/data/arno/labels/exploded_labels/label22.vtk'
    >>> faces, u1,u2, points, u3,u4,u5,u6 = read_vtk(vtk_file)
    >>> order = 3
    >>> zernike_moments(points, faces, order)
    [3535.909590779093,
     82400836.92399396,
     596590.4072532767,
     17284477612.837368,
     66948839.98544357,
     6529913805.079259]

    """
    import numpy as np
    from mindboggle.shapes.zernike.multiproc import MultiprocPipeline

    # Convert lists to numpy arrays:
    if isinstance(points, list):
        points = np.array(points)
    if isinstance(faces, list):
        faces = np.array(faces)

    # Multiprocessor pipeline:
    pl = MultiprocPipeline()

    n_V = 3
    n_faces = pl.size(faces, 0)

    # Geometric moments:
    G = pl.geometric_moments_orig(points, faces, order, n_faces, n_V)

    #
    Z = pl.zernike(G, order)

    # Extract Zernike descriptors:
    descriptors = pl.feature_extraction(Z, order).tolist()

    return descriptors


def zernike_moments_of_largest(points, faces, order=20, exclude_labels=[-1],
                               areas=None):
    """
    Compute the Zernike moments on largest connected segment.

    In case a surface patch is fragmented, we select the largest fragment,
    remove extraneous triangular faces, and reindex indices.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    order : integer
        number of moments to compute
    exclude_labels : list of integers
        labels to be excluded
    areas : numpy array or list of floats (or None)
        surface area scalar values for all vertices

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors of largest connected segment

    Examples
    --------
    >>> # Zernike moments for one label (artificial composite), two fragments:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_of_largest
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> order = 3
    >>> exclude_labels = [-1]
    >>> faces, lines, indices, points, u1, labels, u2,u3 = read_vtk(label_file,
    >>>      return_first=True, return_array=True)
    >>> I19 = [i for i,x in enumerate(labels) if x==19] # pars orbitalis
    >>> I22 = [i for i,x in enumerate(labels) if x==22] # postcentral
    >>> I19.extend(I22)
    >>> faces = remove_faces(faces, I19)
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> #
    >>> zernike_moments_of_largest(points, faces, order, exclude_labels, areas)
    [3535.909590779093,
     82400836.92399396,
     596590.4072532767,
     17284477612.837368,
     66948839.98544357,
     6529913805.079259]
    >>> # View:
    >>> from mindboggle.utils.plots import plot_vtk
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I19] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars')
    >>> plot_vtk(vtk_file)

    """
    import numpy as np

    from mindboggle.utils.segment import select_largest
    from mindboggle.shapes.zernike.zernike import zernike_moments

    if isinstance(areas, list):
        areas = np.array(areas)

    # Check to see if there are enough points:
    min_npoints = order
    npoints = len(points)
    if npoints < min_npoints or len(faces) < min_npoints:
        print("The input size {0} ({1} faces) should be larger "
              "than order {2}".
              format(npoints, len(faces), order))
        return None
    else:

        #---------------------------------------------------------------------
        # Select the largest segment (connected set of indices):
        #---------------------------------------------------------------------
        points, faces = select_largest(points, faces, exclude_labels, areas,
                                       reindex=True)

        # Alert if the number of indices is small:
        if len(points) < min_npoints:
            print("The input size {0} is too small.".format(len(points)))
            return None
        elif faces:

            #-----------------------------------------------------------------
            # Compute Zernike moments:
            #-----------------------------------------------------------------
            descriptors = zernike_moments(points, faces, order)

            return descriptors
        else:
            return None


def zernike_moments_per_label(vtk_file, order=20, exclude_labels=[-1],
                              area_file='', largest_segment=True):
    """
    Compute the Zernike moments per labeled region in a file.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file containing index scalars (labels)
    order : integer
        number of moments to compute
    exclude_labels : list of integers
        labels to be excluded
    area_file :  string
        name of VTK file with surface area scalar values
    largest_segment :  Boolean
        compute moments only for largest segment with a given label?

    Returns
    -------
    descriptors_lists : list of lists of floats
        Zernike descriptors per label
    label_list : list of integers
        list of unique labels for which moments are computed

    Examples
    --------
    >>> # Moments for label 22 (postcentral) in Twins-2-1:
    >>> import os
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> order = 3
    >>> exclude_labels = [0]
    >>> largest_segment = True
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, area_file,
    >>>                           largest_segment)
    ([[3535.909590779093,
       82400836.92399396,
       596590.4072532767,
       17284477612.837368,
       66948839.98544357,
       6529913805.079259]],
     [22])

    """
    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.mesh import remove_faces
    from mindboggle.shapes.zernike.zernike import zernike_moments, \
        zernike_moments_of_largest

    # Read VTK surface mesh file:
    faces, u1,u2, points, u4, labels, u5,u6 = read_vtk(vtk_file)

    # Area file:
    if area_file:
        areas, u1 = read_scalars(area_file)
    else:
        areas = None

    # Loop through labeled regions:
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    descriptors_lists = []
    for label in ulabels:
      if label==22:
        print("DEBUG: COMPUTE FOR ONLY ONE LABEL")

        # Determine the indices per label:
        label_indices = [i for i,x in enumerate(labels) if x == label]
        print('{0} vertices for label {1}'.format(len(label_indices), label))

        # Remove background faces:
        select_faces = remove_faces(faces, label_indices)

        # Compute Zernike moments for the label:
        if largest_segment:
            exclude_labels_inner = [-1]
            descriptors = zernike_moments_of_largest(points, select_faces, order,
                                                     exclude_labels_inner,
                                                     areas)
        else:
            descriptors = zernike_moments(points, select_faces, order)

        # Append to a list of lists of spectra:
        descriptors_lists.append(descriptors)
        label_list.append(label)

    return descriptors_lists, label_list
