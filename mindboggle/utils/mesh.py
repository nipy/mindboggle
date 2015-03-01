#!/usr/bin/env python
"""
Operations on surface mesh vertices.

Authors:
<<<<<<< HEAD
    - Forrest Bao, 2012  (forrest.bao@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
=======
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Bao, 2012  (forrest.bao@gmail.com)
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

<<<<<<< HEAD

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
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> #
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> #
    >>> # Write results to vtk file and view:
    >>> index = 0
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[index] = 1
    >>> IDs[neighbor_lists[index]] = 2
    >>> rewrite_scalars(vtk_file, 'find_neighbors_from_file.vtk', IDs, 'neighbors', IDs)
    >>> plot_surfaces('find_neighbors_from_file.vtk')

    """
    from mindboggle.utils.io_vtk import read_faces_points
    from mindboggle.utils.mesh import find_neighbors

    faces, points, npoints = read_faces_points(input_vtk)

    neighbor_lists = find_neighbors(faces, npoints)

    return neighbor_lists


=======
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
# Find all neighbors from faces
#------------------------------------------------------------------------------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
def find_neighbors(faces, npoints):
    """
    Generate the list of unique, sorted indices of neighboring vertices
    for all vertices in the faces of a triangular mesh.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
<<<<<<< HEAD
=======

>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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
<<<<<<< HEAD
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> faces, points, npoints = read_faces_points(vtk_file)
=======
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, points, npoints = read_faces_points(depth_file)
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    >>> #
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> # Write results to vtk file and view:
    >>> index = 0
    >>> IDs = -1 * np.ones(npoints)
    >>> IDs[index] = 1
    >>> IDs[neighbor_lists[index]] = 2
<<<<<<< HEAD
    >>> rewrite_scalars(vtk_file, 'find_neighbors.vtk', IDs, 'neighbors', IDs)
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> plot_surfaces('find_neighbors.vtk')

    """

    neighbor_lists = [[] for x in range(npoints)]
=======
    >>> rewrite_scalars(depth_file, 'test_find_neighbors.vtk', IDs, 'neighbors', IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_find_neighbors.vtk')

    """

    neighbor_lists = [[] for x in xrange(npoints)]
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

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

<<<<<<< HEAD

=======
#------------------------------------------------------------------------------
# Find neighbors for a given vertex
#------------------------------------------------------------------------------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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

<<<<<<< HEAD

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
    >>> # Extract endpoints from a track in a fold:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file, find_endpoints
    >>> from mindboggle.utils.paths import track_values
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> # Select a single fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> fold_number = 11
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_number]
    >>> # Create a track from the minimum-depth vertex:
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> values, name = read_scalars(vtk_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> seed = indices_fold[np.argmin(values[indices_fold])]
    >>> indices = track_values(seed, indices_fold, neighbor_lists, values, sink=[])
    >>> #
    >>> # Extract endpoints:
    >>> indices_endpoints = find_endpoints(indices, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(len(values))
    >>> IDs[indices_fold] = 1
    >>> IDs[indices] = 2
    >>> IDs[indices_endpoints] = 3
    >>> rewrite_scalars(vtk_file, 'find_endpoints.vtk',
    >>>                 IDs, 'endpoints', IDs)
    >>> plot_surfaces('find_endpoints.vtk')

    """

    # Find vertices with only one neighbor in a set of given indices:
    I = set(indices)
    indices_endpoints = [x for x in indices
                         if len(I.intersection(neighbor_lists[x])) == 1]

    return indices_endpoints


=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
def find_edges(faces):
    """
    Find all edges on a mesh

<<<<<<< HEAD
    Parameters
=======
   Parameters
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
<<<<<<< HEAD
    -------
=======
    --------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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

<<<<<<< HEAD

def find_faces_at_edges(faces):
    """
    For each edge on the mesh, find the two faces that share the edge.

    Parameters
=======
#-----------------------------------------------------------------------------
# find all triangle faces sharing each edge
#-----------------------------------------------------------------------------
def find_faces_at_edges(faces):
    """
    For each edges on the mesh, find the two faces that share the edge.

   Parameters
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
<<<<<<< HEAD
    -------
=======
    --------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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

<<<<<<< HEAD

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
    >>> from mindboggle.utils.mesh import find_faces_with_vertex
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
    >>> from mindboggle.utils.mesh import find_faces_at_vertices
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
    >>> from mindboggle.utils.mesh import find_adjacent_faces
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> find_adjacent_faces(faces)
        [[[-1, 1, 3], [-1, 3, 4]],
         [[-1, 2, 0], [-1, 4, 1]],
         [[4, 3, 1], [1, 1, 2]],
         [[4, 2, 0], [3, 3, 2]],
         [[-1, 3, 2], [-1, 0, 0]]]

    """

    print "Calculating face neighbor list"

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
    >>> from mindboggle.utils.mesh import find_complete_faces
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
=======
#------------------------------------------------------------------------------
# Filter faces
#------------------------------------------------------------------------------
def remove_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices"
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
<<<<<<< HEAD
    indices : integers
        indices to vertices of the surface mesh that are to be retained
=======
    indices : vertex indices to mesh
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

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

<<<<<<< HEAD

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
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_faces_points
    >>> from mindboggle.utils.mesh import reindex_faces_points
    >>> # Reindex faces:
    >>> faces = [[8,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> reindex_faces_points(faces, points=[])
        ([[5, 0, 1], [0, 1, 4], [2, 4, 5], [1, 0, 3]], None)
    >>> # Reindex faces of a single fold of the brain:
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> faces, points, npoints = read_faces_points(fold_file)
    >>> new_faces, new_points, original_indices = reindex_faces_points(faces, points)

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
    >>> from mindboggle.utils.mesh import remove_neighbor_lists
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
    >>> from mindboggle.utils.mesh import reindex_faces_0to1
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
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces, decimate
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> faces, u1,u2, points, u3, scalars, u4,u5 = read_vtk(label_file)
    >>> I22 = [i for i,x in enumerate(labels) if x==14] # postcentral
    >>> faces = remove_faces(faces, I22)
    >>> order = 3
    >>> scale_input = True
    >>> reduction = 0.75
    >>> smooth_steps = 25  # results closer to no decimation than 20 or 30
    >>> save_vtk = False
    >>> output_vtk = ''
    >>> points, faces, scalars, output_vtk = decimate(points, faces, reduction,
    >>>                                               smooth_steps, scalars,
    >>>                                               save_vtk, output_vtk)
    >>> len(points) == 1060
    True
    >>> len(points)
    1060
    >>> # View:
    >>> output_vtk = 'decimated.vtk'
    >>> write_vtk(output_vtk, points, [], [], faces, scalars) # doctest: +SKIP
    >>> plot_surfaces(output_vtk) # doctest: +SKIP

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
    decimate.SetInput(polydata)
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
        smoother.SetInput(decimate.GetOutput())
        smoother.SetNumberOfIterations(smooth_steps)
        smoother.Update()
        out = smoother.GetOutput()
        if save_vtk:
            exporter.SetInput(smoother.GetOutput())
    else:
        decimate.Update()
        out = decimate.GetOutput()
        if save_vtk:
            exporter.SetInput(decimate.GetOutput())

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
                  save_vtk=False, output_vtk=''):
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
    >>> import os
    >>> from mindboggle.utils.mesh import decimate_file
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_vtk = os.path.join(path, 'arno', 'labels', 'label22.vtk')
    >>> #input_vtk='/drop/MB/data/arno/labels/lh.labels.DKT31.manual.vtk'
    >>> save_vtk = True
    >>> output_vtk = ''
    >>> reduction = 0.5
    >>> smooth_steps = 0
    >>> decimate_file(input_vtk, reduction, smooth_steps, save_vtk, output_vtk)
    >>> # View:
    >>> plot_surfaces('decimated.vtk') # doctest: +SKIP

    """
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import decimate

    # Read VTK surface mesh file:
    faces, u1, u2, points, u4, scalars, u5, u6 = read_vtk(input_vtk)

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
    >>> import os
    >>> from mindboggle.utils.mesh import rescale_by_neighborhood
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> indices = []
    >>> nedges = 10
    >>> p = 99
    >>> set_max_to_1 = True
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> background_value = -1
    >>> #
    >>> rescaled_scalars, rescaled_scalars_file = rescale_by_neighborhood(input_vtk,
    >>>     indices, nedges, p, set_max_to_1, save_file, output_filestring, background_value)
    >>> #
    >>> # View rescaled scalar values per fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> rewrite_scalars(rescaled_scalars_file, rescaled_scalars_file,
    >>>                 rescaled_scalars, 'rescaled_depths', folds)
    >>> plot_surfaces(rescaled_scalars_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file, find_neighborhood

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    if not indices:
        indices = [i for i,x in enumerate(scalars) if x != background_value]
    print("  Rescaling {0} scalar values by neighborhood...".format(len(indices)))
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
    >>> import os
    >>> from mindboggle.utils.mesh import rescale_by_label
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> labels_or_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> #
    >>> rescaled_scalars, rescaled_scalars_file = rescale_by_label(input_vtk,
    >>>     labels_or_file, save_file, output_filestring)
    >>> #
    >>> # View rescaled scalar values per fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> rewrite_scalars(rescaled_scalars_file, rescaled_scalars_file,
    >>>                 rescaled_scalars, 'rescaled_depths', folds)
    >>> plot_surfaces(rescaled_scalars_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    print("  Rescaling scalar values within each label...")

    # Load label numbers:
    if isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file, True, True)
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    unique_labels = np.unique(labels)
    unique_labels = [x for x in unique_labels if x >= 0]

    # Loop through labels:
    for label in unique_labels:
        #print("  Rescaling scalar values within label {0} of {1} labels...".format(
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

        rescaled_scalars_file = os.path.join(os.getcwd(), output_filestring + '.vtk')
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
=======
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

def fill_holes(regions, neighbor_lists, exclude_values=[], values=[]):
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
    exclude_values : list of integers
        hole is not filled if any of these values are within the hole
        (prevents cases where surface connected by folds mistaken for holes)
    values : list of integers
        values for vertices, for use in determining which holes to remove

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
    >>> regions = fill_holes(regions, neighbor_lists, [100], values)
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
            # if hole does not include any of the exclude_values
            if exclude_values:
                Ihole = np.where(hole > -1)[0]
                if not len(frozenset(values[Ihole]).intersection(exclude_values)):
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
