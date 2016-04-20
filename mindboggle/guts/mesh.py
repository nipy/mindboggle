#!/usr/bin/env python
"""
Operations on surface mesh vertices.

Authors:
    - Forrest Bao, 2012  (forrest.bao@gmail.com)
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


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
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_mean_curvature'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> neighbor_lists[0:3]
    [[1, 4, 48, 49], [0, 4, 5, 49, 2], [1, 5, 6, 49, 50, 54]]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> index = 100 # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(neighbor_lists)) # doctest: +SKIP
    >>> IDs[index] = 1 # doctest: +SKIP
    >>> IDs[neighbor_lists[index]] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'find_neighbors_from_file.vtk', IDs,
    ...                 'neighbors', IDs) # doctest: +SKIP
    >>> plot_surfaces('find_neighbors_from_file.vtk') # doctest: +SKIP

    """
    from mindboggle.mio.vtks import read_faces_points
    from mindboggle.guts.mesh import find_neighbors

    faces, points, npoints = read_faces_points(input_vtk)

    neighbor_lists = find_neighbors(faces, npoints)

    return neighbor_lists


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
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> npoints = 5
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> neighbor_lists
    [[1, 2, 3, 4], [0, 2, 4, 3], [0, 1, 3], [0, 2, 4, 1], [0, 3, 1]]

    Real example:

    >>> import numpy as np
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.mio.vtks import read_faces_points
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_mean_curvature'])
    >>> faces, points, npoints = read_faces_points(vtk_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> neighbor_lists[0:3]
    [[1, 4, 48, 49], [0, 4, 5, 49, 2], [1, 5, 6, 49, 50, 54]]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> index = 100 # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(neighbor_lists)) # doctest: +SKIP
    >>> IDs[index] = 1 # doctest: +SKIP
    >>> IDs[neighbor_lists[index]] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'find_neighbors.vtk', IDs, 'neighbors', IDs) # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces('find_neighbors.vtk') # doctest: +SKIP

    """

    neighbor_lists = [[] for x in range(npoints)]

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
    >>> from mindboggle.guts.mesh import find_neighbors_vertex
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4]]
    >>> index = 1
    >>> neighbor_lists = find_neighbors_vertex(faces, index)
    >>> neighbor_lists
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


def find_neighborhood(neighbor_lists, indices, nedges=1):
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
    nedges : integer
        number of edges to propagate from indices

    Returns
    -------
    neighborhood : list of integers
        indices to vertices in neighborhood

    Examples
    --------
    >>> from mindboggle.guts.mesh import find_neighborhood
    >>> neighbor_lists = [[0,1],[0,2],[1,4,5],[2],[],[0,1,4,5]]
    >>> indices = [1,3,4]
    >>> neighborhood = find_neighborhood(neighbor_lists, indices, 2)
    >>> neighborhood
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
            seed_list = list(set(local_neighbors).difference(completed))

            # Add to neighborhood:
            neighborhood.extend(seed_list)
            completed.extend(seed_list)

            neighborhood = [int(x) for x in neighborhood]

    return neighborhood


def find_endpoints(indices, neighbor_lists):
    """
    Extract endpoints from connected set of vertices.

    Parameters
    ----------
    indices : list of integers
        indices to connected vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    indices_endpoints : list of integers
        indices to endpoints of connected vertices

    Examples
    --------
    >>> # Find endpoints of fundus in a fold:
    >>> from mindboggle.guts.mesh import find_endpoints
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> from mindboggle.mio.vtks import read_scalars
    >>> urls, fetch_data = prep_tests()
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> fundus_file = fetch_data(urls['left_fundi'])
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fundi, name = read_scalars(fundus_file, True, True)
    >>> background_value = -1
    >>> # Limit number of folds to speed up the test:
    >>> limit_folds = True
    >>> if limit_folds:
    ...     fold_numbers = [2]
    ...     i0 = [i for i,x in enumerate(folds) if x not in fold_numbers]
    ...     folds[i0] = background_value
    ...     fundi[i0] = background_value
    ...     indices = [i for i,x in enumerate(fundi) if x != background_value]
    >>> neighbor_lists = find_neighbors_from_file(fundus_file)
    >>> indices_endpoints = find_endpoints(indices, neighbor_lists)
    >>> indices_endpoints[0:5]
    [68445, 68453]

    View endpoints (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> fundi[indices_endpoints] = 50 # doctest: +SKIP
    >>> rewrite_scalars(fundus_file, 'find_endpoints.vtk', fundi,
    ...     'endpoints', folds, background_value) # doctest: +SKIP
    >>> plot_surfaces('find_endpoints.vtk') # doctest: +SKIP

    """

    # Find vertices with only one neighbor in a set of given indices:
    I = set(indices)
    indices_endpoints = [x for x in indices
                         if len(I.intersection(neighbor_lists[x])) == 1]

    return indices_endpoints


def find_edges(faces):
    """
    Find all edges on a mesh

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    -------
    edges : list of lists of integers
        each element is a 2-tuple of vertex ids representing an edge

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.guts.mesh import find_edges
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


def find_faces_at_edges(faces):
    """
    For each edge on the mesh, find the two faces that share the edge.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    -------
    faces_at_edges : dictionary
        keys are tuples of two vertex IDs and values are 2-tuples of face IDs

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.guts.mesh import find_faces_at_edges
    >>> faces=[[0,1,2], [0,1,4], [1,2,3], [0,2,5]]
    >>> faces_at_edges = find_faces_at_edges(faces)
    >>> faces_at_edges[(0,2)]
    [0, 3]
    >>> faces_at_edges[(2,1)]
    [0, 2]

    Notes ::
        The faces are assumed to be triangular.

    """

    faces_at_edges = {}
    for face_id, face in enumerate(faces):
        for edge in [face[0:2], face[1:3], [face[0], face[2]] ]:
            faces_at_edges.setdefault((edge[0], edge[1]), []).append(face_id)
            faces_at_edges.setdefault((edge[1], edge[0]), []).append(face_id) # make it symmetric

    return faces_at_edges


def find_faces_with_vertex(index, faces):
    """
    For a given vertex, find all faces containing this vertex.
    Note: faces do not have to be triangles.

    Parameters
    ----------
    index : integer
        index to a vertex
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    -------
    faces_with_vertex : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.guts.mesh import find_faces_with_vertex
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> index = 3
    >>> find_faces_with_vertex(index, faces)
    [[0, 2, 3], [0, 3, 4], [4, 3, 1]]

    """
    faces_with_vertex = [x for x in faces if index in x]

    return faces_with_vertex


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
    -------
    faces_at_vertices : list of lists of integers
        faces_at_vertices[i] is a list of faces that contain the i-th vertex

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.guts.mesh import find_faces_at_vertices
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> npoints = 5
    >>> find_faces_at_vertices(faces, npoints)
    [[0, 1, 2, 3], [0, 3, 4], [0, 1], [1, 2, 4], [2, 3, 4]]

    """
    faces_at_vertices = [[] for i in range(npoints)]
    for face_id, face in enumerate(faces):
        for vertex in face:
           faces_at_vertices[vertex].append(face_id)

    return faces_at_vertices


def find_adjacent_faces(faces):
    """
    For each face in a list of faces, find adjacent faces.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
        (-1 indicates no result for a given face or vertex)

    Returns
    -------
    adjacent_faces: list of pairs of lists of three integers
        list 1 indexes three faces adjacent to the three face's edges;
        list 2 indexes three vertices opposite the adjacent faces:
        adjacent_faces[i][0] = [face0, face1, face2], neighbors of face i
        (face0 is the neighbor of face i facing vertex0)
        adjacent_faces[i][1] = [vertex0, vertex1, vertex2] for face i
        (vertex0 is the vertex of face0 not in face i)

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.guts.mesh import find_adjacent_faces
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> adjacent_faces = find_adjacent_faces(faces)
    >>> adjacent_faces[0:2]
    [[[-1, 1, 3], [-1, 3, 4]], [[-1, 2, 0], [-1, 4, 1]]]

    """

    #print("Calculating face neighbor list")

    n_faces = len(faces)

    adjacent_faces = []
    [adjacent_faces.append([[-1,-1,-1], [-1,-1,-1]]) for i in range(n_faces)]

    Done =[]
    [Done.append(0) for i in range(n_faces)]

    # Loop through faces:
    for i1, face1 in enumerate(faces):
        # Loop through remaining faces:
        for i2 in range(i1+1, n_faces):
            face2 = faces[i2]

            # Loop through first two vertices of face:
            for ivertex in [0,1]:
                index1 = face1[ivertex]
                # Loop through remaining vertices of face:
                for index2 in face1[ivertex+1:3]:

                    # If pair of vertices in face2:
                    if index1 in face2 and index2 in face2:

                        # Determine if it is face0, face1 or face2:
                        NbrID1 = 3 - face1.index(index1) - face1.index(index2)
                        NbrID2 = 3 - face2.index(index1) - face2.index(index2)

                        adjacent_faces[i1][0][NbrID1] = i2
                        adjacent_faces[i2][0][NbrID2] = i1
                        adjacent_faces[i1][1][NbrID1] = face2[NbrID2]
                        adjacent_faces[i2][1][NbrID2] = face1[NbrID1]

                        Done[i1] += 1
                        Done[i2] += 1

            # Break if all three neighbors of face1 have been found:
            if Done[i1] == 3:
                break

    return adjacent_faces


def find_complete_faces(indices, faces):
    """
    Given a set of vertices, find the ones that make complete faces.

    Parameters
    ----------
    indices : list of integers
        indices to connected vertices
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    -------
    indices_complete : list of integers
        indices to vertices making up complete faces

    Examples
    --------
    >>> from mindboggle.guts.mesh import find_complete_faces
    >>> faces = [[0,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [3,7,2,5,9,4]
    >>> find_complete_faces(indices, faces)
    [2, 3, 7, 5]

    """

    indices_complete_list = []
    for face in faces:
        if len(list(frozenset(face).intersection(indices))) == 3:
            indices_complete_list.extend(face)
    indices_complete = []
    [indices_complete.append(x) for x in indices_complete_list
     if x not in indices_complete]

    return indices_complete


def keep_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices".

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    indices : list of integers
        indices to vertices of the surface mesh that are to be retained

    Returns
    -------
    faces : list of lists of three integers
        reduced number of faces

    Examples
    --------
    >>> from mindboggle.guts.mesh import keep_faces
    >>> faces = [[1,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [0,1,2,3,4,5]
    >>> keep_faces(faces, indices)
    [[1, 2, 3], [3, 2, 5]]

    """
    import numpy as np

    fs = frozenset(indices)
    faces = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces = np.reshape(np.ravel(faces), (-1, 3))

    #len_faces = len(faces)
    #if verbose and len(faces) < len_faces:
    #    print('Reduced {0} to {1} triangular faces'.
    #        format(len_faces, len(faces)))

    return faces.tolist()


def reindex_faces_points(faces, points=[]):
    """
    Renumber indices in faces and remove points (coordinates) not in faces.

    Parameters
    ----------
    faces : list of lists of integers
        each sublist contains 3 indices of vertices that form a face
        on a surface mesh
    points : list of lists of floats (optional)
        each sublist contains 3-D coordinates of a vertex on a surface mesh

    Returns
    -------
    new_faces : list of lists of integers
        each sublist contains 3 (renumbered) indices of vertices
        that form a face on a surface mesh
    new_points : list of lists of floats
        each (new) sublist contains 3-D coordinates of a vertex on a surface mesh
    original_indices : list integers
        list of indices to original points

    Examples
    --------
    >>> from mindboggle.guts.mesh import reindex_faces_points
    >>> # Reindex faces:
    >>> faces = [[8,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> new_faces, new_points, original_indices = reindex_faces_points(faces,
    ...     points=[])
    >>> new_faces
    [[5, 0, 1], [0, 1, 4], [2, 4, 5], [1, 0, 3]]

    Reindex faces of a limited number of folds of the brain:

    >>> import numpy as np
    >>> from mindboggle.guts.mesh import keep_faces
    >>> from mindboggle.mio.vtks import read_faces_points
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_numbers = [4]
    >>> indices = [i for i,x in enumerate(folds) if x in fold_numbers]
    >>> i0 = [i for i,x in enumerate(folds) if x not in fold_numbers]
    >>> background_value = -1
    >>> folds[i0] = background_value
    >>> faces, points, npoints = read_faces_points(folds_file)
    >>> faces = keep_faces(faces, indices)
    >>> faces[0:3]
    [[51535, 50324, 51529], [50317, 50325, 50326], [50324, 50332, 50333]]
    >>> new_faces, new_points, original_indices = reindex_faces_points(faces,
    ...     points)
    >>> new_faces[0:3]
    [[277, 690, 276], [689, 691, 692], [690, 698, 699]]
    >>> print(np.array_str(np.array(points[0]),
    ...       precision=5, suppress_small=True))
    [-13.7924  -76.0973   -2.57594]
    >>> print(np.array_str(np.array(new_points[0]),
    ...       precision=5, suppress_small=True))
    [-13.7802 -12.3814  57.4042]

    View reindexed fold on surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces
    >>> plot_surfaces('reindex_faces_points.vtk') # doctest: +SKIP

    """
    import numpy as np
    import itertools

    if isinstance(points, list):
        pass
    elif isinstance(points, np.ndarray):
        points = points.tolist()
    else:
        raise IOError("points should be either a list or a numpy array.")

    # set() to remove repeated indices and list() to order them for later use:
    indices_to_keep = list(set(itertools.chain(*faces)))
    reindex = dict([(old_index, new_index)
                    for new_index, old_index in enumerate(indices_to_keep)])

    new_faces = [[reindex[old_index] for old_index in face] for face in faces]

    if points:
        new_points = [points[new_index] for new_index in indices_to_keep]
    else:
        new_points = None

    original_indices = indices_to_keep

    return new_faces, new_points, original_indices


def remove_neighbor_lists(neighbor_lists, indices):
    """
    Remove all but a given set of indices from surface mesh neighbor lists.

    Note :: SLOW!

    Parameters
    ----------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    indices : list of integers
        indices to vertices of the surface mesh

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list has indices to remaining neighboring vertices for each vertex

    Examples
    --------
    >>> from mindboggle.guts.mesh import remove_neighbor_lists
    >>> neighbor_lists = [[1,2,3], [2,3,7], [12,43], [4,7,8], [3,2,5]]
    >>> indices = [0,1,2,3,4,5]
    >>> remove_neighbor_lists(neighbor_lists, indices)
    [[1, 2, 3], [2, 3], [], [4], [2, 3, 5]]

    """

    neighbor_lists = [list(frozenset(indices).intersection(x))
                      for x in neighbor_lists]

    return neighbor_lists


def reindex_faces_0to1(faces):
    """
    Convert 0-indices (Python) to 1-indices (Matlab) for all face indices.

    Parameters
    ----------
    faces : list of lists of integers
        each sublist contains 3 0-indices of vertices that form a face
        on a surface mesh

    Returns
    -------
    faces : list of lists of integers
        each sublist contains 3 1-indices of vertices that form a face
        on a surface mesh

    Examples
    --------
    >>> from mindboggle.guts.mesh import reindex_faces_0to1
    >>> faces = [[0,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> reindex_faces_0to1(faces)
    [[1, 3, 4], [3, 4, 8], [5, 8, 9], [4, 3, 6]]

    """

    faces = [[old_index+1 for old_index in face] for face in faces]

    return faces


def decimate(points, faces, reduction=0.75, smooth_steps=25,
             scalars=[], save_vtk=False, output_vtk=''):
    """
    Decimate vtk triangular mesh with vtk.vtkDecimatePro.

    Parameters
    ----------
    points : list of lists of floats
        each element is a list of 3-D coordinates of a vertex on a surface mesh
    faces : list of lists of integers
        each element is list of 3 indices of vertices that form a face
        on a surface mesh
    reduction : float
        fraction of mesh faces to remove
    smooth_steps : integer
        number of smoothing steps
    scalars : list of integers or floats
        optional scalars for output VTK file
    save_vtk : bool
        output decimated vtk file?
    output_vtk : string
        output decimated vtk file name

    Returns
    -------
    points : list of lists of floats
        decimated points
    faces : list of lists of integers
        decimated faces
    scalars : list of integers or floats
        scalars for output VTK file
    output_vtk : string
        output decimated vtk file

    Examples
    --------
    >>> # Example: Twins-2-1 left postcentral pial surface, 0.75 decimation:
    >>> from mindboggle.guts.mesh import decimate
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_vtk = fetch_data(urls['left_freesurfer_labels'])
    >>> points, f1, f2, faces, scalars, f3, f4, f5 = read_vtk(input_vtk)
    >>> reduction = 0.5
    >>> smooth_steps = 25
    >>> save_vtk = True
    >>> output_vtk = 'decimate.vtk'
    >>> points2, faces2, scalars, output_vtk = decimate(points, faces,
    ...     reduction, smooth_steps, scalars, save_vtk, output_vtk)
    >>> (len(points), len(points2))
    (145069, 72535)
    >>> (len(faces), len(faces2))
    (290134, 145066)

    View decimated surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces('decimate.vtk') # doctest: +SKIP

    """
    import os
    import vtk

    #-------------------------------------------------------------------------
    # vtk points:
    #-------------------------------------------------------------------------
    vtk_points = vtk.vtkPoints()
    [vtk_points.InsertPoint(i, x[0], x[1], x[2]) for i,x in enumerate(points)]

    #-------------------------------------------------------------------------
    # vtk faces:
    #-------------------------------------------------------------------------
    vtk_faces = vtk.vtkCellArray()
    for face in faces:
        vtk_face = vtk.vtkPolygon()
        vtk_face.GetPointIds().SetNumberOfIds(3)
        vtk_face.GetPointIds().SetId(0, face[0])
        vtk_face.GetPointIds().SetId(1, face[1])
        vtk_face.GetPointIds().SetId(2, face[2])
        vtk_faces.InsertNextCell(vtk_face)

    #-------------------------------------------------------------------------
    # vtk scalars:
    #-------------------------------------------------------------------------
    if scalars:
        vtk_scalars = vtk.vtkFloatArray()
        vtk_scalars.SetName("scalars")
        for scalar in scalars:
            vtk_scalars.InsertNextValue(scalar)

    #-------------------------------------------------------------------------
    # vtkPolyData:
    #-------------------------------------------------------------------------
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetPolys(vtk_faces)
    if scalars:
        polydata.GetPointData().SetScalars(vtk_scalars)

    #-------------------------------------------------------------------------
    # Decimate:
    #-------------------------------------------------------------------------
    # We want to preserve topology (not let any cracks form).
    # This may limit the total reduction possible.
    decimate = vtk.vtkDecimatePro()

    # Migrate to VTK6:
    # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
    # Old: decimate.SetInput(polydata)
    decimate.SetInputData(polydata)

    decimate.SetTargetReduction(reduction)
    decimate.PreserveTopologyOn()

    #-------------------------------------------------------------------------
    # Smooth:
    #-------------------------------------------------------------------------
    if save_vtk:
        if not output_vtk:
            output_vtk = os.path.join(os.getcwd(), 'decimated.vtk')
        exporter = vtk.vtkPolyDataWriter()
    else:
        output_vtk = None
    if smooth_steps > 0:
        smoother = vtk.vtkSmoothPolyDataFilter()

        # Migrate to VTK6:
        # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
        # Old: smoother.SetInput(decimate.GetOutput())
        smoother.SetInputConnection(decimate.GetOutputPort())

        smoother.SetNumberOfIterations(smooth_steps)
        smoother.Update()
        out = smoother.GetOutput()

        # Migrate to VTK6:
        # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
        # Old: exporter.SetInput(smoother.GetOutput())
        exporter.SetInputConnection(smoother.GetOutputPort())

    else:
        decimate.Update()
        out = decimate.GetOutput()
        if save_vtk:
            # Migrate to VTK6:
            # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
            # http://stackoverflow.com/questions/29020740/
            #        what-is-the-difference-in-setinputconnection-and-setinput
            # Old: exporter.SetInput(decimate.GetOutput())
            exporter.SetInputConnection(decimate.GetOutputPort())

    #-------------------------------------------------------------------------
    # Export output:
    #-------------------------------------------------------------------------
    if save_vtk:
        exporter.SetFileName(output_vtk)
        exporter.Write()
        if not os.path.exists(output_vtk):
            raise IOError(output_vtk + " not found")

    #-------------------------------------------------------------------------
    # Extract decimated points, faces, and scalars:
    #-------------------------------------------------------------------------
    points = [list(out.GetPoint(point_id))
              for point_id in range(out.GetNumberOfPoints())]
    if out.GetNumberOfPolys() > 0:
        polys = out.GetPolys()
        pt_data = out.GetPointData()
        faces = [[int(polys.GetData().GetValue(j))
                  for j in range(i*4 + 1, i*4 + 4)]
                  for i in range(polys.GetNumberOfCells())]
        if scalars:
            scalars = [pt_data.GetScalars().GetValue(i)
                       for i in range(len(points))]
    else:
        faces = []
        scalars = []

    return points, faces, scalars, output_vtk


def decimate_file(input_vtk, reduction=0.5, smooth_steps=100,
                  save_vtk=True, output_vtk=''):
    """
    Decimate vtk triangular mesh file with vtk.vtkDecimatePro.

    Parameters
    ----------
    input_vtk : string
        input vtk file with triangular surface mesh
    reduction : float
        fraction of mesh faces to remove
    do_smooth : bool
        smooth after decimation?
    save_vtk : bool
        output decimated vtk file?
    output_vtk : string
        output decimated vtk file name

    Returns
    -------
    output_vtk : string
        output decimated vtk file

    Examples
    --------
    >>> from mindboggle.guts.mesh import decimate_file
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_vtk = fetch_data(urls['left_freesurfer_labels'])
    >>> save_vtk = True
    >>> output_vtk = 'decimate.vtk'
    >>> reduction = 0.5
    >>> smooth_steps = 25
    >>> output_vtk = decimate_file(input_vtk, reduction, smooth_steps,
    ...     save_vtk, output_vtk)
    >>> f1, f2, f3, faces1, f4, f5, npoints1, f6 = read_vtk(input_vtk)
    >>> f1, f2, f3, faces2, f4, f5, npoints2, f6 = read_vtk('decimated.vtk')
    >>> (npoints1, npoints2)
    (145069, 72535)
    >>> (len(faces1), len(faces2))
    (290134, 145066)

    View decimated surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces
    >>> plot_surfaces('decimate.vtk') # doctest: +SKIP

    """
    from mindboggle.mio.vtks import read_vtk
    from mindboggle.guts.mesh import decimate

    if not save_vtk:
        raise NotImplementedError()

    # Read VTK surface mesh file:
    points, indices, lines, faces, scalars, scalar_names, npoints, \
            input_vtk = read_vtk(input_vtk)

    # Decimate vtk triangular mesh with vtk.vtkDecimatePro
    points, faces, scalars, output_vtk = decimate(points, faces, reduction,
                                                  smooth_steps, scalars,
                                                  save_vtk, output_vtk)
    return output_vtk


def rescale_by_neighborhood(input_vtk, indices=[], nedges=10, p=99,
    set_max_to_1=True, save_file=False, output_filestring='rescaled_scalars',
    background_value=-1):
    """
    Rescale the scalar values of a VTK file by a percentile value
    in each vertex's surface mesh neighborhood.

    Parameters
    ----------
    input_vtk : string
        name of VTK file with a scalar value for each vertex
    indices : list of integers (optional)
        indices of scalars to normalize
    nedges : integer
        number or edges from vertex, defining the size of its neighborhood
    p : float in range of [0,100]
        percentile used to normalize each scalar
    set_max_to_1 : bool
        set all rescaled values greater than 1 to 1.0?
    save_file : bool
        save output VTK file?
    output_filestring : string (if save_file)
        name of output file
    background_value : integer
        background value

    Returns
    -------
    rescaled_scalars : list of floats
        rescaled scalar values
    rescaled_scalars_file : string (if save_file)
        name of output VTK file with rescaled scalar values

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import rescale_by_neighborhood
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_vtk = fetch_data(urls['left_travel_depth'])
    >>> indices = []
    >>> nedges = 10
    >>> p = 99
    >>> set_max_to_1 = True
    >>> save_file = True
    >>> output_filestring = 'rescale_by_neighborhood'
    >>> background_value = -1
    >>> rescaled, rescaled_file = rescale_by_neighborhood(input_vtk,
    ...     indices, nedges, p, set_max_to_1, save_file, output_filestring,
    ...     background_value)
    >>> scalars1, name = read_scalars(input_vtk)
    >>> print('{0:0.5f}, {1:0.5f}'.format(max(scalars1), max(rescaled)))
    34.95560, 1.00000
    >>> print('{0:0.5f}, {1:0.5f}'.format(np.mean(scalars1), np.mean(rescaled)))
    7.43822, 0.44950

    View rescaled scalar values on surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces(rescaled_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    from mindboggle.guts.mesh import find_neighbors_from_file, find_neighborhood

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    if not indices:
        indices = [i for i,x in enumerate(scalars) if x != background_value]
    #print("  Rescaling {0} scalar values by neighborhood...".format(len(indices)))
    neighbor_lists = find_neighbors_from_file(input_vtk)

    # Loop through vertices:
    rescaled_scalars = scalars.copy()
    for index in indices:

        # Determine the scalars in the vertex's neighborhood:
        neighborhood = find_neighborhood(neighbor_lists, [index], nedges)

        # Compute a high neighborhood percentile to normalize vertex's value:
        normalization_factor = np.percentile(scalars[neighborhood], p)
        rescaled_scalar = scalars[index] / normalization_factor
        rescaled_scalars[index] = rescaled_scalar

    # Make any rescaled value greater than 1 equal to 1:
    if set_max_to_1:
        rescaled_scalars[[x for x in indices if rescaled_scalars[x] > 1.0]] = 1

    rescaled_scalars = rescaled_scalars.tolist()

    #-------------------------------------------------------------------------
    # Return rescaled scalars and file name
    #-------------------------------------------------------------------------
    if save_file:

        rescaled_scalars_file = os.path.join(os.getcwd(), output_filestring + '.vtk')
        rewrite_scalars(input_vtk, rescaled_scalars_file,
                        rescaled_scalars, 'rescaled_scalars', [],
                        background_value)
        if not os.path.exists(rescaled_scalars_file):
            raise IOError(rescaled_scalars_file + " not found")

    else:
        rescaled_scalars_file = None

    return rescaled_scalars, rescaled_scalars_file


def rescale_by_label(input_vtk, labels_or_file, save_file=False,
                     output_filestring='rescaled_scalars',
                     background_value=-1, verbose=False):
    """
    Rescale scalars for each label (such as depth values within each fold).

    Default is to normalize the scalar values of a VTK file by
    a percentile value in each vertex's surface mesh for each label.

    Parameters
    ----------
    input_vtk : string
        name of VTK file with a scalar value for each vertex
    labels_or_file : list or string
        label number for each vertex or name of VTK file with index scalars
    save_file : bool
        save output VTK file?
    output_filestring : string (if save_file)
        name of output file
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    rescaled_scalars : list of floats
        scalar values rescaled for each label, for label numbers not equal to -1
    rescaled_scalars_file : string (if save_file)
        name of output VTK file with rescaled scalar values for each label

    Examples
    --------
    >>> # Rescale depths by neighborhood within each label:
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import rescale_by_label
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_vtk = fetch_data(urls['left_travel_depth'])
    >>> labels_or_file = fetch_data(urls['left_folds'])
    >>> save_file = True
    >>> output_filestring = 'rescale_by_label'
    >>> background_value = -1
    >>> verbose = False
    >>> rescaled, rescaled_label_file = rescale_by_label(input_vtk,
    ...     labels_or_file, save_file, output_filestring, background_value, verbose)
    >>> scalars1, name = read_scalars(input_vtk)
    >>> print('{0:0.5f}, {1:0.5f}'.format(max(scalars1), max(rescaled)))
    34.95560, 1.00000
    >>> print('{0:0.5f}, {1:0.5f}'.format(np.mean(scalars1), np.mean(rescaled)))
    7.43822, 0.30677

    View rescaled scalar values on surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces(rescaled_label_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    from mindboggle.mio.vtks import read_scalars, rewrite_scalars

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    if verbose:
        print("  Rescaling scalar values within each label...")

    # Load label numbers:
    if isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file, True, True)
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    unique_labels = np.unique(labels)

    # Loop through labels:
    for label in unique_labels:
        if verbose:
            print("  Rescaling values within label {0} of {1} labels...".
                format(int(label), len(unique_labels)))
        indices = [i for i,x in enumerate(labels) if x == label]
        if indices:

            # Rescale by the maximum label scalar value:
            scalars[indices] = scalars[indices] / np.max(scalars[indices])
            #print(max(scalars), max(scalars[indices]))

    rescaled_scalars = scalars.tolist()

    #-------------------------------------------------------------------------
    # Return rescaled scalars and file name
    #-------------------------------------------------------------------------
    if save_file:

        rescaled_scalars_file = os.path.join(os.getcwd(),
                                             output_filestring + '.vtk')
        rewrite_scalars(input_vtk, rescaled_scalars_file,
                        rescaled_scalars, 'rescaled_scalars', labels,
                        background_value)
        if not os.path.exists(rescaled_scalars_file):
            raise IOError(rescaled_scalars_file + " not found")

    else:
        rescaled_scalars_file = None

    return rescaled_scalars, rescaled_scalars_file


def area_of_faces(points, faces):
    """
    Compute the areas of all triangles on the mesh.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh

    Returns
    -------
    area: 1-D numpy array
        area[i] is the area of the i-th triangle

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import area_of_faces
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_vtk = fetch_data(urls['left_area'])
    >>> points, f1, f2, faces, f3, f4, f5, f6 = read_vtk(input_vtk)
    >>> area = area_of_faces(points, faces)
    >>> print(np.array_str(area[0:5], precision=5, suppress_small=True))
    [ 0.21703  0.27139  0.29033  0.1717   0.36011]

    """
    import numpy as np

    area = np.zeros(len(faces))

    points = np.array(points)

    for i, triangle in enumerate(faces):

        a = np.linalg.norm(points[triangle[0]] - points[triangle[1]])
        b = np.linalg.norm(points[triangle[1]] - points[triangle[2]])
        c = np.linalg.norm(points[triangle[2]] - points[triangle[0]])
        s = (a+b+c) / 2.0

        area[i] = np.sqrt(s*(s-a)*(s-b)*(s-c))

    return area


def dilate(indices, nedges, neighbor_lists):
    """
    Dilate region on a surface mesh.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to dilate
    nedges : integer
        number of edges to dilate across
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    dilated_indices : list of integers
        indices of original vertices with dilated vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import dilate, find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_travel_depth'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 4
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> dilated_indices = dilate(indices, nedges, neighbor_lists)
    >>> (len(indices), len(dilated_indices))
    (1151, 1545)
    >>> dilated_indices[0:10]
    [50317, 50324, 50325, 50326, 50327, 50332, 50333, 50334, 50339, 50340]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[dilated_indices] = 2 # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'dilate.vtk', IDs, 'dilated_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('dilate.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

    N = find_neighborhood(neighbor_lists, indices, nedges)

    dilated_indices = indices[:]
    dilated_indices.extend(N)

    return dilated_indices


def erode(indices, nedges, neighbor_lists):
    """
    Erode region on a surface mesh.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to erode
    nedges : integer
        number of edges to erode across
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    eroded_indices : list of integers
        indices of original vertices without eroded vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import erode, find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> nedges = 3
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 4
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> eroded_indices = erode(indices, nedges, neighbor_lists)
    >>> (len(indices), len(eroded_indices))
    (1151, 809)

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> IDs[eroded_indices] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'erode.vtk', IDs, 'eroded_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('erode.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

    N1 = find_neighborhood(neighbor_lists, indices, nedges=1)
    N2 = find_neighborhood(neighbor_lists, N1, nedges)

    eroded_indices = list(frozenset(indices).difference(N2))

    return eroded_indices


def extract_edge(indices, neighbor_lists):
    """
    Erode region on a surface mesh to extract the region's edge.

    Parameters
    ----------
    indices : list of integers
        indices of vertices to erode
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    edge_indices : list of integers
        indices of eroded vertices

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import extract_edge
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> # Select a single fold:
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 4
    >>> indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>> edge_indices = extract_edge(indices, neighbor_lists)
    >>> (len(indices), len(edge_indices))
    (1151, 111)

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> IDs = -1 * np.ones(len(folds)) # doctest: +SKIP
    >>> IDs[indices] = 1 # doctest: +SKIP
    >>> IDs[edge_indices] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'extract_edge.vtk', IDs, 'edges_of_fold', IDs) # doctest: +SKIP
    >>> plot_surfaces('extract_edge.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.mesh import find_neighborhood

    N1 = find_neighborhood(neighbor_lists, indices, nedges=1)
    N2 = find_neighborhood(neighbor_lists, N1, nedges=1)
    edge_indices = list(set(N2).intersection(indices))

    return edge_indices


def topo_test(index, values, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    "Simple" is not to be mistaken with the following usage:
    "A vertex is usually assigned one of five possible classifications:
    simple, complex, boundary, interior edge, or corner vertex.
    A simple vertex is surrounded by a closed fan of triangles".

    Parameters
    ----------
    index : integer
        index of vertex
    values : numpy array of integers or floats
        values for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sp : bool
        simple point or not?
    n_inside : integer
        number of neighboring vertices with a value greater than threshold

    Examples
    --------
    >>> # Square with a center vertex:
    >>> # indices [[0,1,2],[3,4,6],[7,8,9]] = 0 and indices [2,4,6] = 1:
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import topo_test
    >>> values = np.array([0,0,1,0,1,0,1,0,0])
    >>> neighbor_lists = [[1,3],[0,2,3,4],[1,4,5],
    ...                   [0,1,4,6],[1,2,3,5,6,7],[2,4,7,8],
    ...                   [3,4,7],[4,5,6,8],[5,7]]
    >>> sps = []
    >>> for index in range(9):
    ...     sp, n_inside = topo_test(index, values, neighbor_lists)
    ...     sps.append(sp)
    >>> sps
    [False, True, True, True, False, True, True, True, False]

    """
    import numpy as np

    # Make sure argument is a numpy array:
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater or less than a class boundary threshold equal to 0.5
    # ("inside" and "outside"); count inside and outside neighbors:
    I_neighbors = neighbor_lists[index]
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > 0.5]
    n_inside = len(inside)
    n_outside = len(I_neighbors) - n_inside

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point:
    if n_outside * n_inside == 0:
        sp = False
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point:
    elif n_outside == 1 or n_inside == 1:
        sp = True
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point:
    else:
        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' indices in a sublist of "N":
        labels = list(range(1, n_inside + 1))
        N = []
        for i_in in range(n_inside):
            new_neighbors = neighbor_lists[inside[i_in]]
            new_neighbors = [x for x in new_neighbors
                             if values[x] > 0.5 if x != index]
            new_neighbors.extend([inside[i_in]])
            N.append(new_neighbors)

        # Consolidate labels of connected vertices --
        # Loop through neighbors (lists within "N"),
        # reassigning the labels for the lists until each label's
        # list(s) has a unique set of vertices:
        change = True
        while change:
            change = False

            # Loop through pairs of inside neighbors
            # and continue if their two labels are different:
            for i in range(n_inside - 1):
                for j in range(i + 1, n_inside):
                    if labels[i] != labels[j]:
                        # Assign the two subsets the same label
                        # if they share at least one vertex,
                        # and continue looping:
                        if set(N[i]).intersection(N[j]):
                            labels[i] = max([labels[i], labels[j]])
                            labels[j] = labels[i]
                            change = True

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other (one unique label):
        D = []
        if len([D.append(x) for x in labels if x not in D]) == 1:
            sp = True
        else:
            sp = False

    return sp, n_inside


# def fill_holes(regions, neighbor_lists, values=[], exclude_range=[],
#                background_value=-1):
#     """
#     Fill holes in regions on a surface mesh by using region boundaries.
#
#     NOTE: assumes one set of connected vertices per region
#
#     Steps ::
#
#         1. Segment region vertex neighbors into connected vertices (region boundaries).
#         2. Remove the largest region boundary, presumably the
#            outer contour of the region, leaving smaller boundaries,
#            presumably the contours of holes within the region.
#         3. Call label_holes() to fill holes with surrounding region numbers.
#
#     Parameters
#     ----------
#     regions : numpy array of integers
#         region numbers for all vertices
#     neighbor_lists : list of lists of integers
#         each list contains indices to neighboring vertices for each vertex
#     values : list of integers
#         values for vertices, for use in determining which holes to remove
#     exclude_range : list of two floats
#         hole is not filled if it contains values within this range
#         (prevents cases where surface connected by folds mistaken for holes)
#     background_value : integer
#         background value
#
#     Returns
#     -------
#     regions : numpy array of integers
#         region numbers for all vertices
#
#     Examples
#     --------
#     >>> import numpy as np
#     >>> from mindboggle.guts.mesh import fill_holes
#     >>> from mindboggle.guts.mesh import find_neighbors_from_file
#     >>> from mindboggle.mio.vtks import read_scalars
#     >>> from mindboggle.mio.fetch_data import prep_tests
#     >>> urls, fetch_data = prep_tests()
#     >>> folds_file = fetch_data(urls['left_folds'])
#     >>> background_value = -1
#     >>> # Select one fold
#     >>> folds, name = read_scalars(folds_file, True, True)
#     >>> fold_number = 4
#     >>> folds[folds != fold_number] = background_value
#     >>> I = np.where(folds==fold_number)[0]
#     >>> neighbor_lists = find_neighbors_from_file(folds_file)
#     >>> ## Find vertex whose removal (with its neighbors) would create a hole:
#     >>> #for index in I:
#     ... #    N1 = neighbor_lists[index]
#     ... #    stop = True
#     ... #    for n in N1:
#     ... #        if any(folds[neighbor_lists[n]] == background_value):
#     ... #            stop = False
#     ... #            break
#     ... #        else:
#     ... #            for f in neighbor_lists[n]:
#     ... #                if any(folds[neighbor_lists[f]] == background_value):
#     ... #                    stop = False
#     ... #                    break
#     ... #    if stop:
#     ... #        break
#     >>> index = I[100]
#     >>> N = neighbor_lists[index]
#     >>> N.append(index)
#     >>> N
#     [36768, 37670, 36769, 37679, 38522, 38529, 37688, 37689, 37678]
#     >>> folds[N] = background_value
#     >>> I = [x for x in I if x not in N]
#
#     View hole (skip test):
#
#     >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
#     >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
#     >>> rewrite_scalars(folds_file, 'hole.vtk', folds, 'hole', folds) # doctest: +SKIP
#     >>> plot_surfaces('hole.vtk') # doctest: +SKIP
#
#     Fill hole:
#
#     >>> exclude_range = []
#     >>> regions = np.copy(folds)
#     >>> values = np.copy(folds)
#     >>> regions = fill_holes(regions, neighbor_lists, values, exclude_range,
#     ...                      background_value)
#     >>> indices = [i for i,x in enumerate(regions) if x != background_value]
#     >>> indices[0:10]
#     [34148, 34149, 34150, 34151, 34152, 34153, 34154, 34155, 34157, 34158]
#
#     View filled hole (skip test):
#
#     >>> rewrite_scalars(folds_file, 'fill_hole.vtk', regions, 'fill_hole', regions) # doctest: +SKIP
#     >>> plot_surfaces('fill_hole.vtk') # doctest: +SKIP
#
#     """
#     import numpy as np
#     from mindboggle.guts.segment import segment
#
#     # Make sure argument is a numpy array
#     if not isinstance(regions, np.ndarray):
#         regions = np.array(regions)
#
#     def label_holes(holes, regions, neighbor_lists):
#         """
#         Fill holes in regions on a surface mesh.
#
#         Parameters
#         ----------
#         holes : list or array of integers
#             hole numbers for all vertices
#         regions : numpy array of integers
#             region numbers for all vertices
#         neighbor_lists : list of lists of integers
#             each list contains indices to neighboring vertices for each vertex
#
#         Returns
#         -------
#         regions : numpy array of integers
#             region numbers for all vertices
#
#         """
#         import numpy as np
#
#         # Make sure argument is a numpy array
#         if not isinstance(regions, np.ndarray):
#             regions = np.array(regions)
#
#         # Identify the vertices for each hole
#         hole_numbers = [x for x in np.unique(holes) if x != background_value]
#         for n_hole in hole_numbers:
#             I = [i for i,x in enumerate(holes) if x == n_hole]
#
#             # Identify neighbors to these vertices
#             N=[]; [N.extend(neighbor_lists[i]) for i in I]
#             if N:
#
#                 # Assign the hole the maximum region ID number of its neighbors
#                 regions[I] = max([regions[x] for x in N])
#
#         return regions
#
#     #-------------------------------------------------------------------------
#     # Find boundaries to holes
#     #-------------------------------------------------------------------------
#     hole_boundaries = background_value * np.ones(len(regions))
#
#     # Identify vertices for each region
#     region_numbers = [x for x in np.unique(regions) if x != background_value]
#     count = 0
#     for n_region in region_numbers:
#         region_indices = np.where(regions == n_region)[0]
#
#         # Identify neighbors to these vertices and their neighbors
#         N = []
#         [N.extend(neighbor_lists[x]) for x in region_indices]
#         N = list(frozenset(N).difference(region_indices))
#         N2 = []
#         [N2.extend(neighbor_lists[x]) for x in N]
#         N.extend(N2)
#         N = list(frozenset(N).difference(region_indices))
#         if N:
#
#             # Segment neighbors into connected vertices (region boundaries)
#             boundaries = segment(N, neighbor_lists)
#
#             # Remove the largest region boundary, presumably the
#             # outer contour of the region, leaving smaller boundaries,
#             # presumably the contours of holes within the region
#             boundary_numbers = [x for x in np.unique(boundaries)
#                                 if x != background_value]
#             max_size = 0
#             max_number = 0
#             for n_boundary in boundary_numbers:
#                 border_indices = np.where(boundaries == n_boundary)[0]
#                 if len(border_indices) > max_size:
#                     max_size = len(border_indices)
#                     max_number = n_boundary
#             boundaries[boundaries == max_number] = background_value
#             boundary_numbers = [x for x in boundary_numbers if x != max_number]
#
#             # Add remaining boundaries to holes array
#             for n_boundary in boundary_numbers:
#                 indices = [i for i,x in enumerate(boundaries) if x == n_boundary]
#                 hole_boundaries[indices] = count
#                 count += 1
#
#     #-------------------------------------------------------------------------
#     # Fill holes
#     #-------------------------------------------------------------------------
#     # If there are any holes
#     if count > 0:
#         hole_numbers = [x for x in np.unique(hole_boundaries)
#                         if x != background_value]
#         background = [i for i,x in enumerate(regions)
#                       if x == background_value]
#
#         # Grow seeds from hole boundaries to fill holes
#         for n_hole in hole_numbers:
#             seed_list = np.where(hole_boundaries == n_hole)[0].tolist()
#             seed_lists = [list(frozenset(background).intersection(seed_list))]
#             hole = segment(background, neighbor_lists, 1, seed_lists)
#
#             # Label the vertices for each hole by surrounding region number
#             # if hole does not include values within exclude_range:
#             if len(exclude_range) == 2:
#                 Ihole = np.where(hole != background_value)[0]
#                 #if not len(frozenset(values[Ihole]).intersection(exclude_range)):
#                 if not [x for x in values[Ihole]
#                         if x > exclude_range[0] if x < exclude_range[1]]:
#                     regions = label_holes(hole, regions, neighbor_lists)
#             else:
#                 regions = label_holes(hole, regions, neighbor_lists)
#
#     return regions


# def close_surface_pair(faces, points1, points2, scalars, background_value=-1):
#     """
#     Close a surface patch by connecting its border vertices with
#     corresponding vertices in a second surface file.
#
#     Assumes no lines or indices when reading VTK files.
#
#     Note ::
#
#         Scalar values different than background define the surface patch.
#         The two sets of points have a 1-to-1 mapping; they are from
#         two surfaces whose corresponding vertices are shifted in position.
#         For pial vs. gray-white matter, the two surfaces are not parallel,
#         so connecting the vertices leads to intersecting faces.
#
#     Parameters
#     ----------
#     faces : list of lists of integers
#         each sublist contains 3 indices of vertices that form a face
#         on a surface mesh
#     points1 : list of lists of floats
#         each sublist contains 3-D coordinates of a vertex on a surface mesh
#     points2 : list of lists of floats
#         points from second surface with 1-to-1 correspondence with points1
#     scalars : numpy array of integers
#         labels used to find foreground vertices
#     background_value : integer
#         scalar value for background vertices
#
#     Returns
#     -------
#     closed_faces : list of lists of integers
#         indices of vertices that form a face on the closed surface mesh
#     closed_points : list of lists of floats
#         3-D coordinates from points1 and points2
#     closed_scalars : list of integers
#         scalar values for points1 and points2
#
#     Examples
#     --------
#     >>> # Build a cube by closing two parallel planes:
#     >>> from mindboggle.guts.mesh import close_surface_pair
#     >>> # Build plane:
#     >>> background_value = -1
#     >>> n = 10  # plane edge length
#     >>> points1 = []
#     >>> for x in range(n):
#     ...     for y in range(n):
#     ...          points1.append([x,y,0])
#     >>> points2 = [[x[0],x[1],1] for x in points1]
#     >>> scalars = [background_value for x in range(len(points1))]
#     >>> p = int(n*(n-1)/2 - 1)
#     >>> for i in [p, p+1, p+n, p+n+1]:
#     ...     scalars[i] = 1
#     >>> faces = []
#     >>> for x in range(n-1):
#     ...     for y in range(n-1):
#     ...         faces.append([x+y*n,x+n+y*n,x+n+1+y*n])
#     ...         faces.append([x+y*n,x+1+y*n,x+n+1+y*n])
#     >>> #write_vtk('plane.vtk', points1, [], [], faces, scalars)
#     >>> #plot_surfaces('plane.vtk')
#     >>> closed_faces, closed_points, closed_scalars = close_surface_pair(faces,
#     ...     points1, points2, scalars, background_value)
#     >>> closed_faces[0:4]
#     [[44, 54, 55], [44, 45, 55], [144, 154, 155], [144, 145, 155]]
#
#     View cube (skip test):
#
#     >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
#     >>> from mindboggle.mio.vtks import write_vtk # doctest: +SKIP
#     >>> write_vtk('cube.vtk', closed_points, [],[], closed_faces,
#     ...     closed_scalars, 'int') # doctest: +SKIP
#     >>> plot_surfaces('cube.vtk') # doctest: +SKIP
#
#     """
#     import sys
#     import numpy as np
#
#     from mindboggle.guts.mesh import find_neighbors, keep_faces
#     from mindboggle.guts.segment import extract_borders
#
#     if isinstance(scalars, list):
#         scalars = np.array(scalars)
#
#     N = len(points1)
#     closed_points = points1 + points2
#
#     # Find all vertex neighbors and surface patch border vertices:
#     neighbor_lists = find_neighbors(faces, N)
#     I = np.where(scalars != background_value)[0]
#     scalars[scalars == background_value] = background_value + 1
#     scalars[I] = background_value + 2
#     scalars = scalars.tolist()
#     borders, u1, u2 = extract_borders(list(range(N)), scalars, neighbor_lists)
#     if not len(borders):
#         sys.exit('There are no border vertices!')
#     borders = [x for x in borders if x in I]
#
#     # Reindex copy of faces and combine with original (both zero-index):
#     indices = list(range(N))
#     indices2 = list(range(N, 2 * N))
#     reindex = dict([(index, indices2[i]) for i, index in enumerate(indices)])
#     faces = keep_faces(faces, I)
#     faces2 = [[reindex[i] for i in face] for face in faces]
#     closed_faces = faces + faces2
#
#     # Connect border vertices between surface patches and add new faces:
#     add_faces = []
#     taken_already = []
#     for index in borders:
#         if index not in taken_already:
#             neighbors = list(set(neighbor_lists[index]).intersection(borders))
#             taken_already.append(index)
#             #taken_already.extend([index] + neighbors)
#             for neighbor in neighbors:
#                 add_faces.append([index, index + N, neighbor])
#                 add_faces.append([index + N, neighbor, neighbor + N])
#     closed_faces = closed_faces + add_faces
#
#     closed_scalars = scalars * 2
#
#     return closed_faces, closed_points, closed_scalars


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules