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
    Use mayavi2 to visualize VTK surface mesh data.

    Inputs
    ------
    vtk_file : string
        name of VTK surface mesh file
    """
    import os, sys

    args = ['mayavi2', '-d', vtk_file, '-m Surface &']
    c = ' '.join(args)
    print(c); os.system(c)

#------------------------------------------------------------------------------
# Plot histogram of VTK surface mesh scalar values
#------------------------------------------------------------------------------
def plot_scalar_histogram(vtk_file, nbins=100):
    """
    Plot histogram of VTK surface mesh scalar values.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.mesh import plot_scalar_histogram
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> plot_scalar_histogram(vtk_file, nbins=500)

    """
    import matplotlib.pyplot as plt
    from mindboggle.utils.io_vtk import read_scalars

    # Load values:
    values, name = read_scalars(vtk_file)

    # Histogram:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(values, nbins, normed=1, facecolor='gray', alpha=0.1)
    plt.show()

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
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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
    faces, lines, indices, points, npoints, scalars, name, input_vtk = read_vtk(vtk_file)

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
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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

#------------------------------------------------------------------------------
# Find neighborhood for given vertices
#------------------------------------------------------------------------------
def find_neighborhood(neighbor_lists, indices, nedges):
    """
    Find neighbors in the neighborhood of given surface mesh vertices.

    For indices to surface mesh vertices, find unique indices for
    vertices in the neighborhood of the vertices.

    Parameters
    ----------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    indices : list of integers
        indices of surface vertices

    Returns
    -------
    neighborhood : list integers
        indices to vertices in neighborhood

    Examples
    --------
    >>> from mindboggle.utils.mesh import find_neighborhood
    >>> neighbor_lists = [[0,1],[0,2],[1,4,5],[2],[],[0,1,4,5]]
    >>> indices = [1,3,4]
    >>> find_neighborhood(neighbor_lists, indices, 2)
        [0, 2, 5]

    """

    # Initialize seed list with indices
    neighborhood = []
    seed_list = indices[:]
    completed = seed_list[:]

    # Propagate nedges away from indices:
    for iedge in range(nedges):

        # Find neighbors of seeds:
        if seed_list:
            local_neighbors = []
            [local_neighbors.extend(neighbor_lists[x]) for x in seed_list]

            # Select neighbors that have not been previously selected:
            seed_list = list(frozenset(local_neighbors).difference(completed))

            # Add to neighborhood:
            neighborhood.extend(seed_list)
            completed.extend(seed_list)

    return neighborhood

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
    indices : integers
        indices to vertices of the surface mesh

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

def renumber_faces(faces, indices):
    """
    Renumber the indices to vertices in faces of a surface mesh.

    This program renumbers the faces lists' indices to vertices,
    so that the new indices are in the range of the number of indices:
    range(len(indices)). This assumes that all of the indices are
    represented in the faces, for example, after running remove_faces().

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    indices : integers
        indices to vertices of the surface mesh, and all contained in faces

    Returns
    -------
    faces_renumbered : list of lists of three integers
        faces with renumbered indices

    Examples
    --------
    >>> from mindboggle.utils.mesh import renumber_faces
    >>> faces = [[8,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [2,3,4,5,8,7]
    >>> renumber_faces(faces, indices)
      [[4, 0, 1], [0, 1, 5], [2, 5, 4], [1, 0, 3]]

    """
    import numpy as np

    faces_ravel = np.ravel(faces)
    faces_renumbered = faces_ravel.copy()

    # Loop through indices:
    for new_index, index in enumerate(indices):
        ifaces = np.where(faces_ravel == index)[0]
        faces_renumbered[ifaces] = new_index

    faces_renumbered = np.reshape(faces_renumbered, (-1, 3))

    return faces_renumbered.tolist()

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
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
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
def skeletonize(binary_array, indices_to_keep, neighbor_lists, values=[]):
    """
    Skeletonize a binary numpy array into 1-vertex-thick curves.
    This does not necessarily find a smooth skeleton.
    It just iteratively removes simple points (see topo_test()),
    optionally in order of lowest to highest values.

    Parameters
    ----------
    binary_array : numpy array of integers
        binary values for all vertices
    indices_to_keep : indices to retain
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    values : list of floats
        optionally remove simple points in order of lowest to highest values

    Returns
    -------
    binary_array : numpy array of integers
        skeleton: binary values for all vertices

    Examples
    --------
    >>> # Extract a skeleton from a fold through a couple of points:
    >>> # (Alternative to connecting vertices with connect_points().
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, \
    >>>                                     read_faces_points, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors, skeletonize
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> # Get neighbor_lists, scalars
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> # Select a single fold:
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> # Test with pre-computed endpoints:
    >>> #endpoints_file = os.path.join(path, 'tests', 'connect_points_test1.vtk')
    >>> #endpoints_file = os.path.join(path, 'tests', 'connect_points_test2.vtk')
    >>> #endpoints, name = read_scalars(endpoints_file)
    >>> #indices_endpoints = [i for i,x in enumerate(endpoints) if x > 1]
    >>> endpoints_file = os.path.join(path, 'tests', 'connect_points_test3.vtk')
    >>> endpoints, name = read_scalars(endpoints_file)
    >>> max_endpoints = max(endpoints)
    >>> indices_endpoints = [i for i,x in enumerate(endpoints) if x == max_endpoints]
    >>> #
    >>> skeleton = skeletonize(fold, indices_endpoints, neighbor_lists)
    >>> #
    >>> # Write out vtk file and view:
    >>> indices_skeleton = [i for i,x in enumerate(skeleton) if x > -1]
    >>> skeleton[indices_endpoints] = 2
    >>> rewrite_scalars(fold_file, 'skeletonize.vtk',
    >>>                 skeleton, 'skeleton', skeleton)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('skeletonize.vtk')

    """
    import numpy as np

    # Make sure argument is a numpy array:
    if not isinstance(binary_array, np.ndarray):
        binary_array = np.array(binary_array)

    # Sort indices to remove simple points in order of lowest to highest values:
    indices = np.where(binary_array)[0]
    if values:
        indices = indices[np.argsort(values[indices])]

    # Loop until all vertices are not simple points:
    exist_simple = True
    while exist_simple == True:
        exist_simple = False

        # For each index:
        for index in indices:

            # Do not update certain indices:
            if binary_array[index] and index not in indices_to_keep:

                # Test to see if index is a simple point:
                update, n_in = topo_test(index, binary_array, neighbor_lists)

                # If a simple point, remove and run again:
                if update and n_in > 1:
                    binary_array[index] = 0
                    exist_simple = True

    return binary_array

#------------------------------------------------------------------------------
# Extract endpoints
#------------------------------------------------------------------------------
def extract_skeleton_endpoints(indices_skeleton, neighbor_lists):
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
    >>> from mindboggle.labels.label import extract_borders, extract_skeleton_endpoints
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_pair_lists = sulcus_boundaries()
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> label_indices = [i for i,x in enumerate(labels) if x in label_pair_lists[0][0]]
    >>> indices_boundaries, label_pairs, foo = extract_borders(label_indices,
    >>>                                               labels, neighbor_lists)
    >>> #
    >>> indices_endpoints = extract_skeleton_endpoints(indices_boundaries, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> end_IDs = -1 * np.ones(len(points))
    >>> end_IDs[indices_boundaries] = 1
    >>> end_IDs[indices_endpoints] = 2
    >>> rewrite_scalars(labels_file, 'test_extract_skeleton_endpoints.vtk',
    >>>                 end_IDs, 'endpoints', end_IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_skeleton_endpoints.vtk')

    """

    # Find vertices with only one connected neighbor
    indices_endpoints = []
    for index in indices_skeleton:
        if len([x for x in neighbor_lists[index] if x in indices_skeleton]) == 1:
            indices_endpoints.append(index)

    return indices_endpoints

#------------------------------------------------------------------------------
# Find points with maximal values that are not too close together.
#------------------------------------------------------------------------------
def find_special_points(points, values, min_directions, min_distance, thr):
    """
    Find 'special' points with maximal values that are not too close together.

    Steps ::

        1. Sort values and find values above the threshold.

        2. Initialize special points with the maximum value,
           remove this value, and loop through the remaining high values.

        3. If there are no nearby special points,
           assign the maximum value vertex as a special point.

    Parameters
    ----------
    points : numpy array of floats
        coordinates for all vertices
    values : list (or array) of integers
        values of some kind to maximize over for all vertices (default -1)
    min_directions : numpy array of floats
        minimum directions for all vertices
    min_distance : integer
        minimum distance
    thr : float
        value threshold in [0,1]

    Returns
    -------
    indices_special : list of integers
        subset of surface mesh vertex indices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_vtk, read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_special_points
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> min_curvature_vector_file = os.path.join(path, 'arno', 'shapes',
    >>>                                          'lh.pial.curv.min.dir.txt')
    >>> faces, lines, indices, points, npoints, values, name, input_vtk = read_vtk(likelihood_file,
    >>>     return_first=True, return_array=True)
    >>> # Select a single fold
    >>> plot_single_fold = True
    >>> if plot_single_fold:
    >>>   fold_ID = 11
    >>>   indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>>   indices_not_fold = [i for i,x in enumerate(folds) if x != fold_ID]
    >>>   values[indices_not_fold] = 0
    >>>   fold_array = -1 * np.ones(len(folds))
    >>>   fold_array[indices_fold] = 1
    >>>   folds = fold_array.tolist()
    >>> #
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> min_distance = 5
    >>> thr = 0.5
    >>> #
    >>> indices_special = find_special_points(points, values, min_directions, min_distance, thr)
    >>> #
    >>> # Write results to vtk file and view:
    >>> values[indices_special] = np.max(values) + 0.1
    >>> rewrite_scalars(likelihood_file, 'find_special_points.vtk',
    >>>                 values, 'special_points_on_values_in_folds', folds)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('find_special_points.vtk')

    """
    import numpy as np
    from operator import itemgetter

    # Make sure arguments are numpy arrays:
    if not isinstance(points, np.ndarray):
        points = np.array(points)
    if not isinstance(min_directions, np.ndarray):
        min_directions = np.array(min_directions)

    max_distance = 2 * min_distance

    # Sort values and find indices for values above the threshold:
    L_table = [[i,x] for i,x in enumerate(values)]
    L_table_sort = np.transpose(sorted(L_table, key=itemgetter(1)))[:, ::-1]
    IL = [int(L_table_sort[0,i]) for i,x in enumerate(L_table_sort[1,:])
          if x > thr]

    # Initialize special points list with the index of the maximum value,
    # remove this value, and loop through the remaining high values:
    if IL:
        indices_special = [IL.pop(0)]
        for imax in IL:

            # Determine if there are any special points
            # near to the current maximum value vertex:
            i = 0
            found = 0
            while i < len(indices_special) and found == 0:

                # Compute Euclidean distance between points:
                D = np.linalg.norm(points[indices_special[i]] - points[imax])

                # If distance less than threshold, consider the point found:
                if D < min_distance:
                    found = 1
                # Compute directional distance between points if they are close:
                elif D < max_distance:
                    dirV = np.dot(points[indices_special[i]] - points[imax],
                                  min_directions[indices_special[i]])
                    # If distance less than threshold, consider the point found:
                    if np.linalg.norm(dirV) < min_distance:
                        found = 1

                i += 1

            # If there are no nearby special points,
            # assign the maximum value vertex as a special point:
            if not found:
                indices_special.append(imax)
    else:
        indices_special = []

    return indices_special

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
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
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
