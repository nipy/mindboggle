#!/usr/bin/env python
"""
Operations on surface mesh vertices.

Authors:
    - Forrest Bao, 2012  (forrest.bao@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    >>> url = urls['left_mean_curvature']
    >>> vtk_file = fetch_data(url)
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> neighbor_lists[0:3]
    [[1, 4, 48, 49], [0, 4, 5, 49, 2], [1, 5, 6, 49, 50, 54]]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> index = 0 # doctest: +SKIP
    >>> IDs = -1 * np.ones(npoints) # doctest: +SKIP
    >>> IDs[index] = 1 # doctest: +SKIP
    >>> IDs[neighbor_lists[index]] = 2 # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'find_neighbors_from_file.vtk', IDs, 'neighbors', IDs) # doctest: +SKIP
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
    >>> url = urls['left_mean_curvature']
    >>> vtk_file = fetch_data(url)
    >>> faces, points, npoints = read_faces_points(vtk_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> neighbor_lists[0:3]
    [[1, 4, 48, 49], [0, 4, 5, 49, 2], [1, 5, 6, 49, 50, 54]]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> index = 0 # doctest: +SKIP
    >>> IDs = -1 * np.ones(npoints) # doctest: +SKIP
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
    >>> from mindboggle.guts.mesh import find_endpoints
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url1 = urls['left_folds']
    >>> folds_file = fetch_data(url1)
    >>> indices = range(1000)
    >>> neighbor_lists = find_neighbors_from_file(folds_file)
    >>> indices_endpoints = find_endpoints(indices, neighbor_lists)
    >>> indices_endpoints
    [848, 886, 920, 987]

    More realistic example -- extract endpoints from a track in a fold:

    >>> import numpy as np
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.paths import track_values
    >>> url2 = urls['left_travel_depth']
    >>> values_file = fetch_data(url2)
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_number]
    >>> values, name = read_scalars(values_file, True, True)
    >>> seed = indices_fold[np.argmin(values[indices_fold])]
    >>> sink = []
    >>> track = track_values(seed, indices_fold, neighbor_lists, values, sink)
    >>> track[0:10]
    [38198, 38207, 39702, 39714, 41177, 41193, 41210, 41222, 41234, 41247]
    >>> indices_endpoints = find_endpoints(track, neighbor_lists)
    >>> indices_endpoints
    [38198, 41325]

    View track in fold on surface (skip test):

    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> T = -1 * np.ones(len(values)) # doctest: +SKIP
    >>> T[indices_fold] = 1 # doctest: +SKIP
    >>> T[track] = 2 # doctest: +SKIP
    >>> T[seed] = 3 # doctest: +SKIP
    >>> T[indices_endpoints] = 4 # doctest: +SKIP
    >>> rewrite_scalars(values_file, 'find_endpoints.vtk', T,
    ...     'endpoints', [], -1) # doctest: +SKIP
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
        adjacent_faces[i]: two lists, each of length 3
        adjacent_faces[i][0] = [face0, face1, face2]:
                                face0 is the neighbor of face i facing vertex0
        adjacent_faces[i][1] = [vertex0, vertex1, vertex2], which is face i:
                                vertex0 is the vertex of face0 not in face i

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
            if Done[i] == 3:
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


def remove_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices".

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    indices : integers
        indices to vertices of the surface mesh that are to be retained

    Returns
    -------
    faces : list of lists of three integers
        reduced number of faces

    Examples
    --------
    >>> from mindboggle.guts.mesh import remove_faces
    >>> faces = [[1,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [0,1,2,3,4,5]
    >>> remove_faces(faces, indices)
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
    >>> new_faces, new_points, original_indices =  reindex_faces_points(faces,
    ...     points=[])
    >>> new_faces
    [[5, 0, 1], [0, 1, 4], [2, 4, 5], [1, 0, 3]]

    Reindex faces of a single fold of the brain:

    >>> import numpy as np
    >>> from mindboggle.mio.vtks import read_faces_points
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url = urls['left_folds']
    >>> folds_file = fetch_data(url)
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_number]
    >>> T = -1 * np.ones(len(folds))
    >>> T[indices_fold] = 1
    >>> output_vtk = rewrite_scalars(folds_file, 'fold.vtk', T, 'fold', [],-1)
    >>> faces, points, npoints = read_faces_points('fold.vtk')
    >>> new_faces, new_points, original_indices = reindex_faces_points(faces,
    ...     points)

    View reindexed fold on surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces('fold.vtk') # doctest: +SKIP

    """
    import itertools

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
    indices : integers
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
    save_vtk : Boolean
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
    >>> url = urls['left_freesurfer_labels']
    >>> input_vtk = fetch_data(url)
    >>> points, f1, f2, faces, scalars, f3, f4, f5 = read_vtk(input_vtk)
    >>> reduction = 0.5
    >>> smooth_steps = 25
    >>> save_vtk = True
    >>> output_vtk = ''
    >>> points2, faces2, scalars, output_vtk = decimate(points, faces,
    ...     reduction, smooth_steps, scalars, save_vtk, output_vtk)
    >>> (len(points), len(points2))
    (145069, 72535)
    >>> (len(faces), len(faces2))
    (290134, 145066)

    View decimated surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces('decimated.vtk') # doctest: +SKIP

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
            raise(IOError(output_vtk + " not found"))

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
    do_smooth : Boolean
        smooth after decimation?
    save_vtk : Boolean
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
    >>> url = urls['left_freesurfer_labels']
    >>> input_vtk = fetch_data(url)
    >>> save_vtk = True
    >>> output_vtk = ''
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

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces('decimated.vtk') # doctest: +SKIP

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
    set_max_to_1 : Boolean
        set all rescaled values greater than 1 to 1.0?
    save_file : Boolean
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
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url = urls['left_travel_depth']
    >>> input_vtk = fetch_data(url)
    >>> indices = []
    >>> nedges = 10
    >>> p = 99
    >>> set_max_to_1 = True
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> background_value = -1
    >>> rescaled, rescaled_file = rescale_by_neighborhood(input_vtk,
    ...     indices, nedges, p, set_max_to_1, save_file, output_filestring,
    ...     background_value)
    >>> scalars1, name = read_scalars(input_vtk)
    >>> (max(scalars1), max(rescaled))
    (34.9556, 1.0)
    >>> (np.mean(scalars1), np.mean(rescaled))
    (7.4382172059361906, 0.44950180785631405)

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
                        rescaled_scalars, 'rescaled_scalars')
        if not os.path.exists(rescaled_scalars_file):
            raise(IOError(rescaled_scalars_file + " not found"))

    else:
        rescaled_scalars_file = None

    return rescaled_scalars, rescaled_scalars_file


def rescale_by_label(input_vtk, labels_or_file, save_file=False,
                     output_filestring='rescaled_scalars'):
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
    save_file : Boolean
        save output VTK file?
    output_filestring : string (if save_file)
        name of output file

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
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url1 = urls['left_travel_depth']
    >>> url2 = urls['left_folds']
    >>> input_vtk = fetch_data(url1)
    >>> labels_vtk = fetch_data(url2)
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> rescaled, rescaled_label_file = rescale_by_label(input_vtk,
    ...     labels_vtk, save_file, output_filestring)
    >>> scalars1, name = read_scalars(input_vtk)
    >>> (max(scalars1), max(rescaled))
    (34.9556, 1.0)
    >>> (np.mean(scalars1), np.mean(rescaled))
    (7.4382172059361906, -0.049475798634033216)

    View rescaled scalar values on surface (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces(rescaled_label_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    from mindboggle.mio.vtks import read_scalars, rewrite_scalars

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    #print("  Rescaling scalar values within each label...")

    # Load label numbers:
    if isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file, True, True)
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    unique_labels = np.unique(labels)
    unique_labels = [x for x in unique_labels if x >= 0]

    # Loop through labels:
    for label in unique_labels:
        #print("  Rescaling values within label {0} of {1} labels...".format(
        #    int(label), len(unique_labels)))
        indices = [i for i,x in enumerate(labels) if x == label]
        if indices:

            # Rescale by the maximum label scalar value:
            scalars[indices] = scalars[indices] / np.max(scalars[indices])

    rescaled_scalars = scalars.tolist()

    #-------------------------------------------------------------------------
    # Return rescaled scalars and file name
    #-------------------------------------------------------------------------
    if save_file:

        rescaled_scalars_file = os.path.join(os.getcwd(),
                                             output_filestring + '.vtk')
        rewrite_scalars(input_vtk, rescaled_scalars_file,
                        rescaled_scalars, 'rescaled_scalars', labels)
        if not os.path.exists(rescaled_scalars_file):
            raise(IOError(rescaled_scalars_file + " not found"))

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
    >>> from mindboggle.guts.mesh import area_of_faces
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> url = urls['left_area']
    >>> input_vtk = fetch_data(url)
    >>> points, f1, f2, faces, f3, f4, f5, f6 = read_vtk(input_vtk)
    >>> area = area_of_faces(points, faces)
    >>> area[0:5]
    array([ 0.21703181,  0.27138698,  0.29033419,  0.17170339,  0.36010545])

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


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()