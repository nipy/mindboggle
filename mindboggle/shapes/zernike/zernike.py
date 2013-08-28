#!/usr/bin/python
"""
Compute the Zernike moments of a collection of points.


Authors:
    - Arthur Mikhno, 2013, Columbia University (original MATLAB code)
    - Brian Rossa, 2013, Tank Think Labs, LLC (port to Python)
    - Arno Klein, 2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def zernike_moments(points, faces, order=10, scale_input=True):
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
    scale_input : Boolean
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors

    Examples
    --------
    >>> # Example 1: simple cube:
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> points = [[0,0,0], [1,0,0], [0,0,1], [0,1,1], [1,0,1], [0,1,0], [1,1,1], [1,1,0]]
    >>> faces = [[0,2,4], [0,1,4], [2,3,4], [3,4,5], [3,5,6], [0,1,7]]
    >>> order = 3
    >>> scale_input = True
    >>> zernike_moments(points, faces, order, scale_input)
    [0.0918881492369654,
     0.09357431096617608,
     0.04309029164656885,
     0.06466432586854755,
     0.03820155248327533,
     0.04138011726544602]
    >>> # Example 2: simple cube (with inner diagonal plane):
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> #path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = '/drop/cube.vtk'
    >>> faces, u1,u2, points, u3,u4,u5,u6 = read_vtk(vtk_file)
    >>> order = 3
    >>> scale_input = True
    >>> zernike_moments(points, faces, order, scale_input)
    [0.0,
     1.5444366221695725e-21,
     0.0,
     2.081366518964347e-21,
     5.735003646768394e-05,
     2.433866250546253e-21]
    >>> # Example 3: superior frontal pial surface patch:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, labels, u4,u5 = read_vtk(label_file)
    >>> I28 = [i for i,x in enumerate(labels) if x==28] # superior frontal
    >>> faces = remove_faces(faces, I28)
    >>> order = 3
    >>> scale_input = True
    >>> zernike_moments(points, faces, order, scale_input)
    [0.01566588675857814,
     0.012346130475603009,
     0.018993428582801408,
     0.022470635336483052,
     0.016625931229056885,
     0.014155674468859897]

    """
    import numpy as np

    from mindboggle.utils.mesh import reindex_faces_0to1
    from mindboggle.shapes.zernike.multiproc import MultiprocPipeline

    # Convert 0-indices (Python) to 1-indices (Matlab) for all face indices:
    index1 = True
    if index1:
        faces = reindex_faces_0to1(faces)

    # Convert lists to numpy arrays:
    if isinstance(points, list):
        points = np.array(points)
    if isinstance(faces, list):
        faces = np.array(faces)

    # Translate all points so that they are centered at their mean,
    # and scale them so that they are bounded by a unit sphere:
    if scale_input:
        center = np.mean(points, axis=0)
        points = points - center
        maxd = np.max(np.sqrt(np.sum(points**2, axis=1)))
        points = points / maxd

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


def zernike_moments_of_largest(points, faces, order=10, exclude_labels=[-1],
                               areas=None, scale_input=True):
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
    scale_input : Boolean
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())

    Returns
    -------
    descriptors : list of floats
        Zernike descriptors of largest connected segment

    Examples
    --------
    >>> # Example 1: superior frontal pial surface patch:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_of_largest
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> order = 3
    >>> exclude_labels = [-1]
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> scale_input = True
    >>> faces, lines, indices, points, u1, labels, u2,u3 = read_vtk(label_file)
    >>> I28 = [i for i,x in enumerate(labels) if x==28] # superior frontal
    >>> faces = remove_faces(faces, I28)
    >>> zernike_moments_of_largest(points, faces, order, exclude_labels, areas, scale_input)
    [0.01566588675857814,
     0.012346130475603009,
     0.018993428582801408,
     0.022470635336483052,
     0.016625931229056885,
     0.014155674468859897]
    # Example 2: superior frontal + pars triangularis pial surface patches:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_of_largest
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> order = 3
    >>> exclude_labels = [-1]
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> scale_input = True
    >>> faces, lines, indices, points, u1, labels, u2,u3 = read_vtk(label_file)
    >>> I28 = [i for i,x in enumerate(labels) if x==28] # superior frontal
    >>> I20 = [i for i,x in enumerate(labels) if x==20] # pars triangularis
    >>> I28.extend(I20)
    >>> faces = remove_faces(faces, I28)
    >>> zernike_moments_of_largest(points, faces, order, exclude_labels, areas, scale_input)
    [0.01566588675857814,
     0.012346130475603009,
     0.018993428582801408,
     0.022470635336483052,
     0.016625931229056885,
     0.014155674468859897]
    >>> # View two surface patches:
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I28] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars')
    >>> plot_surfaces(vtk_file)

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
            descriptors = zernike_moments(points, faces, order, scale_input)

            return descriptors
        else:
            return None


def zernike_moments_per_label(vtk_file, order=10, exclude_labels=[-1],
                              area_file='', largest_segment=True,
                              close_file='', do_decimate=False,
                              reduction=0.5, smooth_steps=100,
                              background_value=-1, scale_input=True):
    """
    Compute the Zernike moments per labeled region in a file.

    Optionally close each label's surface by connecting its borders to
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
    background_value : integer
        background value
    scale_input : Boolean
        translate and scale each object so it is bounded by a unit sphere?
        (this is the expected input to zernike_moments())

    Returns
    -------
    descriptors_lists : list of lists of floats
        Zernike descriptors per label
    label_list : list of integers
        list of unique labels for which moments are computed

    Examples
    --------
    >>> # Example 1: superior frontal pial surface patch -- no decimation:
    >>> import os
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> order = 3
    >>> exclude_labels = [0]
    >>> largest_segment = True
    >>> close_file = ''
    >>> do_decimate = False
    >>> reduction = 0.5
    >>> smooth_steps = 100
    >>> background_value = -1
    >>> scale_input = True
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, area_file,
    >>>     largest_segment, close_file, do_decimate, reduction, smooth_steps, background_value, scale_input)
    ([[0.01566588675857814,
       0.012346130475603009,
       0.018993428582801408,
       0.022470635336483052,
       0.016625931229056885,
       0.014155674468859897]],
     [28])
    >>> # Example 2: superior frontal pial surface patch -- WITH decimation:
    >>> import os
    >>> from mindboggle.shapes.zernike.zernike import zernike_moments_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> order = 3
    >>> exclude_labels = [0]
    >>> largest_segment = True
    >>> close_file = ''
    >>> do_decimate = True
    >>> reduction = 0.5
    >>> smooth_steps = 100
    >>> background_value = -1
    >>> scale_input = True
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, area_file,
    >>>     largest_segment, close_file, do_decimate, reduction, smooth_steps, background_value, scale_input)
    ([[0.02453764381289278,
       0.03697607218757975,
       0.00422907681893017,
       0.009325945207818874,
       0.014734760259388037,
       0.005008394221210818]],
     [28])
    >>> # Example 3: superior frontal pial -- surface closing + decimation:
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
    >>> background_value = -1
    >>> scale_input = True
    >>> zernike_moments_per_label(vtk_file, order, exclude_labels, area_file,
    >>>     largest_segment, close_file, do_decimate, reduction, smooth_steps, background_value, scale_input)
    ([[0.03660348461019181,
       0.055324987603587394,
       0.005843199772889861,
       0.013885747292519622,
       0.022316298827496545,
       0.006667325690449897]],
     [28])

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
      #if label==28:
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
            L = background_value * np.ones(npoints)
            L[Ilabel] = 1
            new_faces, new_points, u1 = close_surfaces(faces, points, points2,
                                                       L, background_value)
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
            exclude_labels_inner = [background_value]
            descriptors = zernike_moments_of_largest(new_points, new_faces,
                order, exclude_labels_inner, areas, scale_input)
        else:
            descriptors = zernike_moments(new_points, new_faces, order,
                                          scale_input)

        #---------------------------------------------------------------------
        # Append to a list of lists of spectra:
        #---------------------------------------------------------------------
        descriptors_lists.append(descriptors)
        label_list.append(label)

    return descriptors_lists, label_list
