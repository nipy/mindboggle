#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(points, faces, order=20, index1=True):
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
    index1 : Boolean
        Convert 0-indices (Python) to 1-indices (Matlab) for all face indices?

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
    >>> index1 = True
    >>> zernike_moments(points, faces, order, index1)
    [0.03978873577297383,
     0.07597287228593756,
     0.08322411181893402,
     0.20233902300388348,
     0.117303689366221,
     0.1457966728239256]
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
    >>> index1 = True
    >>> zernike_moments(points, faces, order, index1)
    [7562.751480397972,
     143262239.5171249,
     1107670.7893994227,
     28487908892.820065,
     112922387.17238183,
     10250734140.30357]

    """
    import numpy as np

    from mindboggle.utils.mesh import reindex_faces_0to1
    from mindboggle.shapes.zernike.multiproc import MultiprocPipeline

    # Convert 0-indices to 1-indices for all face indices:
    if index1:
        faces = reindex_faces_0to1(faces)

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
    [7562.751480397972,
     143262239.5171249,
     1107670.7893994227,
     28487908892.820065,
     112922387.17238183,
     10250734140.30357]
    >>> # View two fragments:
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
                              area_file='', largest_segment=True,
                              close_file='', do_decimate=False,
                              reduction=0.5, smooth_steps=100):
    """
    Compute the Zernike moments per labeled region in a file.

    Optionally close each label's surface by connecting to its borders to
    those of a second corresponding surface, and decimate the resulting mesh.

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
    close_file : string
        second corresponding surface if closing (ex: gray-white surface)
    do_decimate : Boolean
        decimate each label's mesh?
    reduction : float
        fraction of mesh faces to remove for decimation
    smooth_steps : integer
        number of smoothing steps for decimation

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
    >>> close_file = os.path.join(path, 'arno', 'freesurfer', 'lh.white.vtk')
    >>> do_decimate = True
    >>> reduction = 0.5
    >>> smooth_steps = 100
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, area_file,
    >>>     largest_segment, close_file, do_decimate, reduction, smooth_steps)

    """
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk, read_scalars, read_points, write_vtk
    from mindboggle.utils.mesh import remove_faces, close_surfaces, decimate
    from mindboggle.shapes.zernike.zernike import zernike_moments, \
        zernike_moments_of_largest

    #-------------------------------------------------------------------------
    # Read VTK surface mesh file and area file:
    #-------------------------------------------------------------------------
    faces, u1,u2, points, npoints, labels, u4,u5 = read_vtk(vtk_file)

    if area_file:
        areas, u1 = read_scalars(area_file)
    else:
        areas = None

    #-------------------------------------------------------------------------
    # Loop through labeled regions:
    #-------------------------------------------------------------------------
    ulabels = []
    [ulabels.append(int(x)) for x in labels if x not in ulabels
     if x not in exclude_labels]
    label_list = []
    descriptors_lists = []
    for label in ulabels:
      #if label==11:
      #  print("DEBUG: COMPUTE FOR ONLY ONE LABEL")

        #---------------------------------------------------------------------
        # Determine the indices per label:
        #---------------------------------------------------------------------
        Ilabel = [i for i,x in enumerate(labels) if x == label]
        print('  {0} vertices for label {1}'.format(len(Ilabel), label))

        #---------------------------------------------------------------------
        # Close surface:
        #---------------------------------------------------------------------
        if close_file:

            # Read second VTK surface mesh file:
            points2 = read_points(close_file)
            L = -1 * np.ones(npoints)
            L[Ilabel] = 1
            new_faces, new_points, u1 = close_surfaces(faces, points, points2,
                                                       L, background_value=-1)
            #write_vtk('c'+str(label)+'.vtk', new_points, [],[], new_faces, u1)

        #---------------------------------------------------------------------
        # Remove background faces:
        #---------------------------------------------------------------------
        else:
            new_faces = remove_faces(faces, Ilabel)
            new_points = points

        #---------------------------------------------------------------------
        # Decimate surface:
        #---------------------------------------------------------------------
        if do_decimate:
            new_points, new_faces, u1,u2 = decimate(new_points, new_faces,
                reduction, smooth_steps, [], save_vtk=False)
            #write_vtk('d'+str(label)+'.vtk', new_points, [],[], new_faces, [])

        #---------------------------------------------------------------------
        # Compute Zernike moments for the label:
        #---------------------------------------------------------------------
        if largest_segment:
            exclude_labels_inner = [-1]
            descriptors = zernike_moments_of_largest(new_points, new_faces,
                order, exclude_labels_inner, areas)
        else:
            descriptors = zernike_moments(new_points, new_faces, order)

        #---------------------------------------------------------------------
        # Append to a list of lists of spectra:
        #---------------------------------------------------------------------
        descriptors_lists.append(descriptors)
        label_list.append(label)

    return descriptors_lists, label_list
