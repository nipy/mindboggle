#!/usr/bin/env python
"""
Operations on surface mesh vertices.

Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Bao, 2012  (forrest.bao@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Plot VTK surface mesh
#------------------------------------------------------------------------------
def plot_vtk(vtk_file):
    """
    Use mayavi2 to visualize .vtk surface mesh data.

    Inputs
    ------
    vtk_file : string
        name of VTK surface mesh file
    """
    import os, sys

    args = ['mayavi2', '-d' , vtk_file, '-m Surface &']
    c = ' '.join(args)
    print(c); os.system(c)

#------------------------------------------------------------------------------
# Apply affine transform to the points of a VTK surface mesh
#------------------------------------------------------------------------------
def read_itk_transform(affine_transform_file):
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
    affine_transform_file : string
        name of ITK affine transform file

    Returns
    -------
    affine_transform : numpy array
        4x4 affine transform matrix
    fixed_parameters : numpy array
        FixedParameters vector

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.mesh import read_itk_transform
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> affine_transform_file = os.path.join(path, 'arno', 'mri',
    >>>                             't1weighted_brain.MNI152Affine.txt')
    >>> read_itk_transform(affine_transform_file)
        (array([[  9.07680e-01,   4.35290e-02,   1.28917e-02, -7.94889e-01],
               [ -4.54455e-02,   8.68937e-01,   4.06098e-01, -1.83346e+01],
               [  1.79439e-02,  -4.30013e-01,   7.83074e-01, -3.14767e+00],
               [  0.00000e+00,   0.00000e+00,   0.00000e+00, 1.00000e+00]]),
         [-0.60936, 21.1593, 10.6148])

    """
    import numpy as np

    affine_transform = np.eye(4)

    # Read ITK transform file
    fid = open(affine_transform_file, 'r')
    affine_lines = fid.readlines()

    transform = affine_lines[3]
    transform = transform.split()
    transform = [np.float(x) for x in transform[1::]]
    transform = np.reshape(transform, (4,3))
    linear_transform = transform[0:3,:]
    translation = transform[3,:]
    affine_transform[0:3,0:3] = linear_transform
    affine_transform[0:3,3] = translation

    fixed_parameters = affine_lines[4]
    fixed_parameters = fixed_parameters.split()
    fixed_parameters = [np.float(x) for x in fixed_parameters[1::]]

    return affine_transform, fixed_parameters

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
    >>> from mindboggle.utils.mesh import apply_affine_transform, plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> transform_file = os.path.join(path, 'arno', 'mri',
    >>>                               't1weighted_brain.MNI152Affine.txt')
    >>> vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> apply_affine_transform(transform_file, vtk_file)
    >>> # View
    >>> plot_vtk('affine_lh.pial.depth.vtk')


    """
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk, write_vtk
    from mindboggle.utils.mesh import read_itk_transform

    print("\n\n\nUNDER CONSTRUCTION!!!\n\n\n")

    # Read ITK affine transform file
#    transform, fixed_parameters = read_itk_transform(transform_file)

    import nibabel as nb
    subject_path = '/Applications/freesurfer/subjects/Twins-2-1'
    native_volume_mgz = subject_path + '/mri/orig/001.mgz'
    conformed_volume_mgz = subject_path + '/mri/brain.mgz'
    M = np.array([[-1,0,0,128],
                  [0,0,1,-128],
                  [0,-1,0,128],
                  [0,0,0,1]],dtype=float)
    native = nb.freesurfer.load(native_volume_mgz)
    conformed = nb.freesurfer.load(conformed_volume_mgz)
    affine_native = native.get_affine()
    affine_conformed = conformed.get_affine()
    transform = np.dot(affine_conformed, np.linalg.inv(M))

    # Read VTK file
    faces, lines, indices, points, npoints, scalars, name = read_vtk(vtk_file)

    # Transform points
    points = np.array(points)
#    points += fixed_parameters

    points = np.concatenate((points, np.ones((np.shape(points)[0],1))), axis=1)
    affine_points = np.transpose(np.dot(transform, np.transpose(points)))[:,0:3]
#    affine_points -= fixed_parameters
#    #affine_points += [0,256,0]

    # Output transformed VTK file
    output_file = os.path.join(os.getcwd(), 'affine_' + os.path.basename(vtk_file))

    # Write VTK file
    write_vtk(output_file, affine_points.tolist(), indices, lines, faces, scalars, name)

    return affine_points, output_file

#------------------------------------------------------------------------------
# Find all neighbors from faces in a VTK mesh file
#------------------------------------------------------------------------------
def find_neighbors_from_file(input_vtk):
    """
    Generate the list of unique, sorted indices of neighboring vertices
    for all vertices in the faces of a triangular mesh in a VTK file.

    Parameters
    ----------
    input_vtk : string
        name of input VTK file containing surface mesh

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.io_vtk import rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> #
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> #
    >>> # Write results to vtk file and view:
    >>> index = 0
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[index] = 1
    >>> IDs[neighbor_lists[index]] = 2
    >>> rewrite_scalars(depth_file, 'test_find_neighbors_from_file.vtk', IDs, 'neighbors', IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_find_neighbors_from_file.vtk')

    """
    from mindboggle.utils.io_vtk import read_faces_points
    from mindboggle.utils.mesh import find_neighbors

    faces, points, npoints = read_faces_points(input_vtk)

    neighbor_lists = find_neighbors(faces, npoints)

    return neighbor_lists

#------------------------------------------------------------------------------
# Find all neighbors from faces
#------------------------------------------------------------------------------
def find_neighbors(faces, npoints):
    """
    Generate the list of unique, sorted indices of neighboring vertices
    for all vertices in the faces of a triangular mesh.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    npoints: integer
        number of vertices on the mesh

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> npoints = 5
    >>> find_neighbors(faces, npoints)
        [[1, 2, 3, 4], [0, 2, 4, 3], [0, 1, 3], [0, 2, 4, 1], [0, 3, 1]]

    >>> # Real example:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.utils.io_vtk import read_faces_points, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> #
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> # Write results to vtk file and view:
    >>> index = 0
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[index] = 1
    >>> IDs[neighbor_lists[index]] = 2
    >>> rewrite_scalars(depth_file, 'test_find_neighbors.vtk', IDs, 'neighbors', IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_find_neighbors.vtk')

    """

    neighbor_lists = [[] for x in xrange(npoints)]

    for face in faces:
        [v0, v1, v2] = face
        if v1 not in neighbor_lists[v0]:
            neighbor_lists[v0].append(v1)
        if v2 not in neighbor_lists[v0]:
            neighbor_lists[v0].append(v2)

        if v0 not in neighbor_lists[v1]:
            neighbor_lists[v1].append(v0)
        if v2 not in neighbor_lists[v1]:
            neighbor_lists[v1].append(v2)

        if v0 not in neighbor_lists[v2]:
            neighbor_lists[v2].append(v0)
        if v1 not in neighbor_lists[v2]:
            neighbor_lists[v2].append(v1)

    return neighbor_lists

#------------------------------------------------------------------------------
# Find neighbors for a given vertex
#------------------------------------------------------------------------------
def find_neighbors_vertex(faces, index):
    """
    Find neighbors to a surface mesh vertex.

    For a set of surface mesh faces and the index of a surface vertex,
    find unique indices for neighboring vertices.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    index : int
        index of surface vertex

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Examples
    --------
    >>> from mindboggle.utils.mesh import find_neighbors_vertex
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4]]
    >>> index = 1
    >>> find_neighbors_vertex(faces, index)
        [0, 2, 4]

    """
    import numpy as np

    # Make sure argument is a numpy array
    if not isinstance(faces, np.ndarray):
        faces = np.array(faces)

    # Create list of vertex indices sharing the same faces as "index"
    I = [faces[np.where(faces[:,i] == index)[0], :] for i in (0,1,2)]

    # Create single list from nested lists
    I = [int(x) for lst in I for sublst in lst for x in sublst]

    # Find unique indices not equal to "index"
    neighbor_list = []; [neighbor_list.append(x)
                         for x in I if x not in neighbor_list if x != index]

    return neighbor_list

#-----------------------------------------------------------------------------
# find all triangle faces centered at each node on the mesh
#-----------------------------------------------------------------------------
def find_faces_at_vertices(faces, npoints):
    """
    For each vertex, find all faces containing this vertex.
    Note: faces do not have to be triangles.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    npoints: integer
        number of vertices on the mesh

    Returns
    --------
    faces_at_vertex : list of lists of integers
        faces_at_vertices[i] is a list of faces that contain the i-th vertex

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh import find_faces_at_vertices
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> npoints = 5
    >>> find_faces_at_vertices(faces, npoints)
        [[0, 1, 2, 3], [0, 3, 4], [0, 1], [1, 2, 4], [2, 3, 4]]

    """
    faces_at_vertices = [[] for i in xrange(npoints)]
    for face_id, face in enumerate(faces):
        for vertex in face:
           faces_at_vertices[vertex].append(face_id)

    return faces_at_vertices

#-----------------------------------------------------------------------------
# find all edges on the mesh
#-----------------------------------------------------------------------------
def find_edges(faces):
    """
    Find all edges on a mesh

   Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    --------
    edges : list of lists of integers
        each element is a 2-tuple of vertex ids representing an edge

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh import find_edges
    >>> faces=[[0,1,2], [0,1,4], [1,2,3], [0,2,5]]
    >>> find_edges(faces)
    [[0, 1], [1, 2], [0, 2], [1, 4], [0, 4], [2, 3], [1, 3], [2, 5], [0, 5]]

    """
    edges = [ ]
    for face in faces:
        for edge in [face[0:2], face[1:3], [face[0], face[2]] ]:
            if not edge in edges: # I know that this is costly
                edges.append(edge)

    return edges

#-----------------------------------------------------------------------------
# find all triangle faces sharing each edge
#-----------------------------------------------------------------------------
def find_faces_at_edges(faces):
    """
    For each edges on the mesh, find the two faces that share the edge.

   Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    --------
    faces_at_edges : dictionary
        keys are tuples of two vertex IDs and values are 2-tuples of face IDs

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh import find_faces_at_edges
    >>> faces=[[0,1,2], [0,1,4], [1,2,3], [0,2,5]]
    >>> find_faces_at_edges(faces)
        {(0, 1): [0, 1],
         (0, 2): [0, 3],
         (0, 4): [1],
         (0, 5): [3],
         (1, 0): [0, 1],
         (1, 2): [0, 2],
         (1, 3): [2],
         (1, 4): [1],
         (2, 0): [0, 3],
         (2, 1): [0, 2],
         (2, 3): [2],
         (2, 5): [3],
         (3, 1): [2],
         (3, 2): [2],
         (4, 0): [1],
         (4, 1): [1],
         (5, 0): [3],
         (5, 2): [3]}

    Notes ::
        The faces are assumed to be triangular.

    """

    faces_at_edges = {}
    for face_id, face in enumerate(faces):
        for edge in [face[0:2], face[1:3], [face[0], face[2]] ]:
            faces_at_edges.setdefault((edge[0], edge[1]), []).append(face_id)
            faces_at_edges.setdefault((edge[1], edge[0]), []).append(face_id) # make it symmetric

    return faces_at_edges

#------------------------------------------------------------------------------
# Filter faces
#------------------------------------------------------------------------------
def remove_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices"

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    indices : vertex indices to mesh

    Returns
    -------
    faces : list of lists of three integers
        reduced number of faces

    Examples
    --------
    >>> from mindboggle.utils.mesh import remove_faces
    >>> faces = [[1,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [0,1,2,3,4,5]
    >>> remove_faces(faces, indices)
      Reduced 4 to 2 triangular faces.
      [[1, 2, 3], [3, 2, 5]]

    """
    import numpy as np

    len_faces = len(faces)
    fs = frozenset(indices)
    faces = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces = np.reshape(np.ravel(faces), (-1, 3))
    if len(faces) < len_faces:
        print('Reduced {0} to {1} triangular faces'.format(len_faces, len(faces)))

    return faces.tolist()

#------------------------------------------------------------------------------
# Fill holes
#------------------------------------------------------------------------------
def label_holes(holes, regions, neighbor_lists):
    """
    Fill holes in regions on a surface mesh.

    Parameters
    ----------
    holes : list or array of integers
        hole numbers for all vertices (default -1)
    regions : numpy array of integers
        region numbers for all vertices (default -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    regions : numpy array of integers
        region numbers for all vertices (default -1)

    """
    import numpy as np

    # Make sure argument is a numpy array
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    # Identify the vertices for each hole
    hole_numbers = [x for x in np.unique(holes) if x > -1]
    for n_hole in hole_numbers:
        I = [i for i,x in enumerate(holes) if x == n_hole]

        # Identify neighbors to these vertices
        N=[]; [N.extend(neighbor_lists[i]) for i in I]
        if N:

            # Assign the hole the maximum region ID number of its neighbors
            regions[I] = max([regions[x] for x in N])

    return regions

def fill_holes(regions, neighbor_lists, values=[], exclude_range=[]):
    """
    Fill holes in regions on a surface mesh by using region boundaries.

    NOTE: assumes one set of connected vertices per region

    Steps ::

        1. Segment region vertex neighbors into connected vertices (region boundaries).
        2. Remove the largest region boundary, presumably the
           outer contour of the region, leaving smaller boundaries,
           presumably the contours of holes within the region.
        3. Call label_holes() to fill holes with surrounding region numbers.

    Parameters
    ----------
    regions : numpy array of integers
        region numbers for all vertices (default -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    values : list of integers
        values for vertices, for use in determining which holes to remove
    exclude_range : list of two floats
        hole is not filled if it contains values within this range
        (prevents cases where surface connected by folds mistaken for holes)

    Returns
    -------
    regions : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors, remove_faces, fill_holes
    >>> from mindboggle.labels.segment import segment
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #
    >>> # Select one fold
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> n_fold = np.unique(folds)[1]
    >>> folds[folds != n_fold] = -1
    >>> #
    >>> # Make two holes in fold (values of -1 and excluded values)
    >>> # Hole 1:
    >>> # Find a vertex whose removal (with its neighbors) would create a hole
    >>> I = np.where(folds==n_fold)[0]
    >>> for index1 in I:
    >>>     N1 = neighbor_lists[index1]
    >>>     stop = True
    >>>     for n in N1:
    >>>         if any(folds[neighbor_lists[n]] == -1):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == -1):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index1] = -1
    >>> folds[N1] = -1
    >>> # Hole 2:
    >>> I = np.where(folds==n_fold)[0]
    >>> for index2 in I:
    >>>     N2 = neighbor_lists[index2]
    >>>     stop = True
    >>>     for n in N2:
    >>>         if any(folds[neighbor_lists[n]] == -1):
    >>>             stop = False
    >>>             break
    >>>         else:
    >>>             for f in neighbor_lists[n]:
    >>>                 if any(folds[neighbor_lists[f]] == -1):
    >>>                     stop = False
    >>>                     break
    >>>     if stop:
    >>>         break
    >>> folds[index2] = -1
    >>> folds[N2] = -1
    >>> values = np.zeros(len(folds))
    >>> values[index2] = 100
    >>> values[N2] = 200
    >>> #
    >>> # Write holes to vtk file and view:
    >>> holes = folds.copy()
    >>> holes[index1] = 10
    >>> holes[N1] = 20
    >>> holes[index2] = 30
    >>> holes[N2] = 40
    >>> indices = [i for i,x in enumerate(holes) if x > -1]
    >>> write_vtk('test_holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), [holes.tolist()], ['holes'])
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_holes.vtk')
    >>> #
    >>> # Fill Hole 1 but not Hole 2:
    >>> # (because values has an excluded value in the hole)
    >>> regions = np.copy(folds)
    >>> regions = fill_holes(regions, neighbor_lists, values, [99,101])
    >>> #
    >>> # Write results to vtk file and view:
    >>> indices = [i for i,x in enumerate(regions) if x > -1]
    >>> write_vtk('test_fill_holes.vtk', points, indices, lines,
    >>>           remove_faces(faces, indices), regions.tolist(), 'regions')
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_fill_holes.vtk')

    """
    import numpy as np
    from mindboggle.utils.mesh import label_holes
    from mindboggle.labels.segment import segment

    # Make sure argument is a numpy array
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    #--------------------------------------------------------------------------
    # Find boundaries to holes
    #--------------------------------------------------------------------------
    hole_boundaries = -1 * np.ones(len(regions))

    # Identify vertices for each region
    region_numbers = [x for x in np.unique(regions) if x > -1]
    count = 0
    for n_region in region_numbers:
        region_indices = np.where(regions == n_region)[0]

        # Identify neighbors to these vertices a]nd their neighbors
        N = []
        [N.extend(neighbor_lists[x]) for x in region_indices]
        N = list(frozenset(N).difference(region_indices))
        N2 = []
        [N2.extend(neighbor_lists[x]) for x in N]
        N.extend(N2)
        N = list(frozenset(N).difference(region_indices))
        if N:

            # Segment the neighbors into connected vertices (region boundaries)
            boundaries = segment(N, neighbor_lists)

            # Remove the largest region boundary, presumably the
            # outer contour of the region, leaving smaller boundaries,
            # presumably the contours of holes within the region
            boundary_numbers = [x for x in np.unique(boundaries) if x > -1]
            max_size = 0
            max_number = 0
            for n_boundary in boundary_numbers:
                boundary_indices = np.where(boundaries == n_boundary)[0]
                if len(boundary_indices) > max_size:
                    max_size = len(boundary_indices)
                    max_number = n_boundary
            boundaries[boundaries == max_number] = -1
            boundary_numbers = [x for x in boundary_numbers if x != max_number]

            # Add remaining boundaries to holes array
            for n_boundary in boundary_numbers:
                indices = [i for i,x in enumerate(boundaries) if x == n_boundary]
                hole_boundaries[indices] = count
                count += 1

    #--------------------------------------------------------------------------
    # Fill holes
    #--------------------------------------------------------------------------
    # If there are any holes
    if count > 0:
        hole_numbers = [x for x in np.unique(hole_boundaries) if x > -1]
        background = [i for i,x in enumerate(regions) if x == -1]

        # Grow seeds from hole boundaries to fill holes (background value -1)
        for n_hole in hole_numbers:
            seed_list = np.where(hole_boundaries == n_hole)[0].tolist()
            seed_lists = [list(frozenset(background).intersection(seed_list))]
            hole = segment(background, neighbor_lists, 1, seed_lists)

            # Label the vertices for each hole by surrounding region number
            # if hole does not include values within exclude_range:
            if len(exclude_range) == 2:
                Ihole = np.where(hole > -1)[0]
                #if not len(frozenset(values[Ihole]).intersection(exclude_range)):
                if not [x for x in values[Ihole]
                        if x > exclude_range[0] and x < exclude_range[1]]:
                    regions = label_holes(hole, regions, neighbor_lists)
            else:
                regions = label_holes(hole, regions, neighbor_lists)

    return regions

#------------------------------------------------------------------------------
# Test for simple points
#------------------------------------------------------------------------------
def topo_test(index, values, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    Parameters
    ----------
    index : index of vertex
    values : numpy array of integers or floats
        values for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sp : simple point or not?: Boolean
    n_inside : number of neighboring vertices greater than threshold

    """
    import numpy as np

    # Make sure argument is a numpy array
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater or less than the class boundary threshold for HMMF (0.5)
    # ("inside" and "outside"); count the number of inside and outside neighbors
    I_neighbors = neighbor_lists[index]
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > 0.5]
    n_inside = len(inside)
    n_outside = len(I_neighbors) - n_inside

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point
    if n_outside * n_inside == 0:
        sp = False
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point
    elif n_outside == 1 or n_inside == 1:
        sp = True
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point
    else:
        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' indices in a sublist of "N"
        labels = range(1, n_inside + 1)
        N = []
        for i_in in range(n_inside):
            new_neighbors = neighbor_lists[inside[i_in]]
            new_neighbors = [x for x in new_neighbors
                             if values[x] > 0.5 if x != index]
            new_neighbors.extend([inside[i_in]])
            N.append(new_neighbors)

        # Consolidate labels of connected vertices:
        # Loop through neighbors (lists within "N"),
        # reassigning the labels for the lists until each label's
        # list(s) has a unique set of vertices
        change = True
        while change:
            change = False

            # Loop through pairs of inside neighbors
            # and continue if their two labels are different
            for i in range(n_inside - 1):
                for j in range(i + 1, n_inside):
                    if labels[i] != labels[j]:
                        # Assign the two subsets the same label
                        # if they share at least one vertex,
                        # and continue looping
                        if frozenset(N[i]).intersection(N[j]):
                            labels[i] = max([labels[i], labels[j]])
                            labels[j] = labels[i]
                            change = True

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other (one unique label)
        D = []
        if len([D.append(x) for x in labels if x not in D]) == 1:
            sp = True
        else:
            sp = False

    return sp, n_inside


#------------------------------------------------------------------------------
# Skeletonize
#------------------------------------------------------------------------------
def skeletonize(binary_array, indices_to_keep, neighbor_lists):
    """
    Skeletonize a binary numpy array into 1-vertex-thick curves.
    This does not necessarily find a smooth skeleton.
    It just iteratively removes simple points (see topo_test()).

    Parameters
    ----------
    binary_array : numpy array of integers
        binary values for all vertices
    indices_to_keep : indices to retain
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    binary_array : numpy array of integers
        skeleton: binary values for all vertices

    Examples
    --------
    >>> # Extract a skeleton from a fold through a couple of points:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors, skeletonize
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> faces, lines, indices, points, npoints, folds, name = read_vtk(folds_file,
    >>>     return_first=True, return_array=True)
    >>> n_fold = max(folds)
    >>> folds[folds != n_fold] = -1
    >>> indices_fold = [i for i,x in enumerate(folds) if x > -1]
    >>> indices = [indices_fold[0], indices_fold[-1]]
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> skeleton = skeletonize(folds, indices, neighbor_lists)
    >>> #
    >>> # Write out vtk file and view:
    >>> rewrite_scalars(folds_file, 'test_skeletonize.vtk',
    >>>                 skeleton, 'skeleton', skeleton)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_skeletonize.vtk')

    """
    import numpy as np

    # Make sure argument is a numpy array
    if not isinstance(binary_array, np.ndarray):
        binary_array = np.array(binary_array)

    # Loop until all vertices are not simple points
    indices = np.where(binary_array)[0]
    exist_simple = True
    while exist_simple == True:
        exist_simple = False

        # For each index
        for index in indices:

            # Do not update certain indices
            if binary_array[index] and index not in indices_to_keep:

                # Test to see if index is a simple point
                update, n_in = topo_test(index, binary_array, neighbor_lists)

                # If a simple point, remove and run again
                if update and n_in > 1:
                    binary_array[index] = 0
                    exist_simple = True

    return binary_array

#------------------------------------------------------------------------------
# Extract endpoints
#------------------------------------------------------------------------------
def extract_endpoints(indices_skeleton, neighbor_lists):
    """
    Extract endpoints from connected set of vertices

    Parameters
    ----------
    indices_skeleton : list of integers
        indices to connected vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    indices_endpoints : list of integers
        indices to endpoints of connected vertices

    Examples
    --------
    >>> # Extract endpoints from a sulcus label boundary segment
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.label import extract_borders, extract_endpoints
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_pair_lists = sulcus_boundaries()
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> label_indices = [i for i,x in enumerate(labels) if x in label_pair_lists[0][0]]
    >>> indices_boundaries, label_pairs, foo = extract_borders(label_indices,
    >>>                                               labels, neighbor_lists)
    >>> #
    >>> indices_endpoints = extract_endpoints(indices_boundaries, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> end_IDs = -1 * np.ones(len(points))
    >>> end_IDs[indices_boundaries] = 1
    >>> end_IDs[indices_endpoints] = 2
    >>> rewrite_scalars(labels_file, 'test_extract_endpoints.vtk',
    >>>                 end_IDs, 'endpoints', end_IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_endpoints.vtk')

    """

    # Find vertices with only one connected neighbor
    indices_endpoints = []
    for index in indices_skeleton:
        if len([x for x in neighbor_lists[index] if x in indices_skeleton]) == 1:
            indices_endpoints.append(index)

    return indices_endpoints

#------------------------------------------------------------------------------
# Example: watershed() segmentation
#------------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    import numpy as np
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.labels.segment import watershed, propagate
    from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
        return_first=True, return_array=True)

    indices = np.where(depths > 0.11)[0]  # high to speed up
    neighbor_lists = find_neighbors(faces, npoints)

    segments = watershed(depths, indices, neighbor_lists,
                         depth_ratio=0.1, tolerance=0.01, remove_fraction=0.75)

    # Write results to vtk file and view:
    rewrite_scalars(depth_file, 'test_segment.vtk',
                    segments, 'segments', segments)
    from mindboggle.utils.mesh import plot_vtk
    plot_vtk('test_segment.vtk')
