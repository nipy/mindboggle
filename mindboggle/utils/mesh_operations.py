#!/usr/bin/python
"""
Find neighbors to a surface mesh vertex.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------
# Find all neighbors from file
#-----------------------------
def find_all_neighbors_from_file(surface_file):
    """
    Generate the list of unique, sorted indices of neighboring vertices
    for all vertices in the faces of a triangular mesh file.

    Parameters
    ----------
    surface_file : str
        name of VTK file

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices

    """
    from utils.io_vtk import load_scalar

    points, faces, scalars = load_scalar(surface_file)

    neighbor_lists = find_all_neighbors(faces)

    return neighbor_lists

#------------------------------
# Find all neighbors from faces
#------------------------------
def find_all_neighbors(faces):
    """
    Generate the list of unique, sorted indices of neighboring vertices
    for all vertices in the faces of a triangular mesh.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices

    Example
    -------
    >>> from extract.identify_sulci import generate_neighbor_lists
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> find_all_neighbors(faces)
        [[1, 2, 3, 4], [0, 2, 4, 3], [1, 0, 3], [2, 0, 4, 1], [0, 3, 1]]
    >>> from utils.io_vtk import load_scalar
    >>> points, faces, scalars = load_scalar('lh.pial.depth.vtk')
    >>> neighbor_lists = find_all_neighbors(faces)

    """
    import numpy as np

    n_vertices = np.max(faces) + 1
    neighbor_lists = [[] for i in xrange(n_vertices)]

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

#----------------------------------
# Find neighbors for a given vertex
#----------------------------------
def find_neighbors(faces, index):
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
    neighbor_list : list of integers
        indices of neighboring vertices

    Example
    -------
    >>> from extract.fundi_hmmf.find_points import find_neighbors
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4]]
    >>> index = 1
    >>> find_neighbors(faces, index)
        [0, 2, 4]

    """

    import numpy as np

    # Make sure argument is a numpy array
    if type(faces) != np.ndarray:
        faces = np.array(faces)

    # Create list of vertex indices sharing the same faces as "index"
    I = [faces[np.where(faces[:,i] == index)[0], :] for i in (0,1,2)]

    # Create single list from nested lists
    I = [int(x) for lst in I for sublst in lst for x in sublst]

    # Find unique indices not equal to "index"
    neighbor_list = []; [neighbor_list.append(x)
                         for x in I if x not in neighbor_list if x != index]

    return neighbor_list

#===================
# Find anchor points
#===================
def find_anchors(vertices, L, min_directions, min_distance, thr):
    """
    Find anchor points.

    Assign maximum likelihood vertices as "anchor points"
    while ensuring that the anchor points are not close to one another.
    Anchor points are used to construct curves.

    Note: only the sulcus end points are strictly necessary
          (so that the fundus doesn't shrink)

    Parameters
    ----------
    vertices : [#vertices x 3] numpy array
    L : fundus likelihood values: [#vertices x 1] numpy array
    min_directions : [#vertices x 1] numpy array
    min_distance : minimum distance
    thr : likelihood threshold

    Returns
    -------
    anchors : list of subset of surface mesh vertex indices

    """

    import numpy as np
    from operator import itemgetter

    # Make sure arguments are numpy arrays
    if type(vertices) != np.ndarray:
        vertices = np.array(vertices)
    if type(L) != np.ndarray:
        L = np.array(L)
    if type(min_directions) != np.ndarray:
        min_directions = np.array(min_directions)

    max_distance = 2 * min_distance

    # Sort likelihood values and find indices for values above the threshold
    L_table = [[i,x] for i,x in enumerate(L)]
    L_table_sort = np.transpose(sorted(L_table, key=itemgetter(1)))[:, ::-1]
    IL = [int(L_table_sort[0,i]) for i,x in enumerate(L_table_sort[1,:])
          if x > thr]

    # Initialize anchors list with the index of the maximum likelihood value,
    # remove this value, and loop through the remaining high likelihoods
    anchors = [IL.pop(0)]
    for imax in IL:

        # Determine if there are any anchor points
        # near to the current maximum likelihood vertex
        i = 0
        found = 0
        while i < len(anchors) and found == 0:

            # Compute Euclidean distance between points
            D = np.linalg.norm(vertices[anchors[i], :] - vertices[imax, :])

            # If distance less than threshold, consider the point found
            if D < min_distance:
                found = 1
            # Compute directional distance between points if they are close
            elif D < max_distance:
                dirV = np.dot(vertices[anchors[i], :] - vertices[imax, :],
                              min_directions[anchors[i], :])
                # If distance less than threshold, consider the point found
                if np.linalg.norm(dirV) < min_distance:
                    found = 1

            i += 1

        # If there are no nearby anchor points,
        # assign the maximum likelihood vertex as an anchor point
        if not found:
            anchors.append(imax)

    return anchors

#========================
# Segment surface patches
#========================
def segment_surface(faces, seeds, neighbor_lists, n_vertices,
                    min_seeds, min_patch_size):
    """
    Segment a surface into contiguous patches (seed region growing).

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    seeds : mesh vertex indices for vertices to be segmented
            list or [#seeds x 1] array
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices (see return_arrays)
    n_vertices : #vertices total (seeds are a subset)
    min_seeds : minimum number of seeds (vertices) per triangle for inclusion
    min_patch_size : minimum size of segmented set of vertices

    Returns
    -------
    segments : segment indices for patches: [#seeds x 1] numpy array
    n_segments : #segments
    max_patch_label : index for largest segmented set of vertices

    """
    import numpy as np

    # Initialize segments and seeds (indices of deep vertices)
    segments = np.zeros(n_vertices)
    n_seeds = len(seeds)

    # Remove faces with fewer than min_seeds seeds to speed up computation
    fs = frozenset(seeds)
    faces_seeds = [lst for lst in faces
                   if len(fs.intersection(lst)) >= min_seeds]
    faces_seeds = np.reshape(np.ravel(faces_seeds), (-1, 3))
    print('    Reduced {0} to {1} faces.'.format(len(faces),
                                                 len(faces_seeds)))

    # Loop until all seed vertices segmented
    print('    Grow {0} seed vertices...'.format(n_seeds))
    max_patch_size = 0
    max_patch_label = 1
    n_segments = 0
    counter = 0
    TEMP0 = np.zeros(n_vertices)
    while n_seeds >= min_patch_size:
        TEMP = np.copy(TEMP0)

        # Select a seed vertex (selection does not affect result)
        #I = [seeds[round(np.random.rand() * (n_seeds - 1))]]
        I = [seeds[0]]
        # Region grown from seed
        Ipatch = I[:]

        # Grow region about the seed vertex until
        # there are no more connected seed vertices available.
        # Continue loop if there are newly selected neighbors.
        loop = 1
        while loop:
            loop = 0
            TEMP[I] = 1
            Inew = []

            # Find neighbors for each selected seed vertex
            for index in I:
                neighbors = neighbor_lists[index]
                # Select neighbors that have not been previously selected
                if len(neighbors) > 0:
                    neighbors = [x for x in neighbors if TEMP[x] == 0]
                    TEMP[neighbors] = 2
                    Inew.extend(neighbors)

                    # Continue looping
                    loop = 1
            I = Inew
            Ipatch.extend(Inew)

        # Disregard vertices already visited
        seeds = list(frozenset(seeds).difference(Ipatch))
        n_seeds = len(seeds)

        # Assign counter number to segmented patch
        # if patch size is greater than min_patch_size
        size_patch = len(Ipatch)
        if size_patch >= min_patch_size:
            counter += 1
            n_segments = counter
            segments[Ipatch] = n_segments

            # Display current number and size of patch
            if size_patch > 1:
                print('    Segmented patch {0}: {1} vertices. {2} seeds remaining...'.
                      format(n_segments, size_patch, n_seeds))

            # Find the maximum patch size
            if size_patch > max_patch_size:
                max_patch_size = size_patch
                max_patch_label = counter

    return segments, n_segments, max_patch_label


#-----------
# Fill holes
#-----------
def fill_holes(faces, patches, holes, n_holes, neighbor_lists):
    """
    Fill holes in surface mesh patches.

    Parameters
    ----------
    faces : surface mesh vertex indices [#faces x 3] numpy array
    patches : [#vertices x 1] numpy array
    holes : [#vertices x 1] numpy array
    n_holes : [#vertices x 1] numpy array
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices (see return_arrays)

    Returns
    -------
    patches : [#vertices x 1] numpy array

    """

    import numpy as np

    # Make sure arguments are numpy arrays
    if type(faces) != np.ndarray:
        faces = np.array(faces)
    if type(patches) != np.ndarray:
        patches = np.array(patches)

    # Look for vertices that have a fold label and are
    # connected to any of the vertices in the current hole,
    # and assign the hole the maximum label number
    for i in range(1, n_holes + 1):
        indices_holes = np.where(holes == i)[0]
        # Loop until a labeled neighbor is found
        for index in indices_holes:
            # Find neighboring vertices to the hole
            neighbors = neighbor_lists[index]
            # If there are any neighboring labels,
            # assign the hole the maximum label
            # of its neighbors and end the while loop
            for fold_neighbor in patches[neighbors]:
                if fold_neighbor > 0:
                    patches[indices_holes] = fold_neighbor
                    break

    return patches

#-----------------------
# Test for simple points
#-----------------------
def simple_test(faces, index, values, thr, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    Parameters
    ----------
    faces : [#faces x 3] numpy array
    index : index of vertex
    values : values: [#vertices x 1] numpy array
    thr : threshold
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    sp : simple point or not?: Boolean
    n_inside : number of neighboring vertices greater than threshold

    """

    import numpy as np

    # Make sure arguments are numpy arrays
    if type(faces) != np.ndarray:
        faces = np.array(faces)
    if type(values) != np.ndarray:
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater than a threshold, thr,
    # and those less than or equal to thr ("inside" and "outside").
    # Count the number of "inside" and "outside" neighbors
    I_neighbors = neighbor_lists[index]
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > thr]
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
                             if values[x] > thr if x != index]
            new_neighbors.extend([inside[i_in]])
            N.append(new_neighbors)

        # Consolidate labels of connected vertices:
        # Loop through neighbors (lists within "N"),
        # reassigning the labels for the lists until each label's
        # list(s) has a unique set of vertices
        change = 1
        while change > 0:
            change = 0

            # Loop through pairs of inside neighbors
            # and continue if their two labels are different
            for i in range(n_inside - 1):
                for j in range(i + 1, n_inside):
                    if labels[i] != labels[j]:
                        # Assign the two subsets the same label
                        # if they share at least one vertex,
                        # and continue looping
                        if len(frozenset(N[i]).intersection(N[j])) > 0:
                            labels[i] = max([labels[i], labels[j]])
                            labels[j] = labels[i]
                            change = 1

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other (one unique label)
        D = []
        if len([D.append(x) for x in labels if x not in D]) == 1:
            sp = True
        else:
            sp = False

    return sp, n_inside


#------------
# Skeletonize
#------------
def skeletonize(B, indices_to_keep, faces, neighbor_lists):
    """
    Skeletonize a binary numpy array into 1-vertex-thick curves.

    Parameters
    ----------
    B : binary [#vertices x 1] numpy array
    indices_to_keep : indices to retain
    faces : indices of triangular mesh vertices: [#faces x 3] numpy array
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    B : binary skeleton: numpy array

    """

    import numpy as np

    # Make sure arguments are numpy arrays
    if type(B) != np.ndarray:
        B = np.array(B)
    if type(faces) != np.ndarray:
        faces = np.array(faces)

    # Loop until all vertices are not simple points
    indices = np.where(B)[0]
    exist_simple = True
    while exist_simple == True:
        exist_simple = False

        # For each index
        for index in indices:

            # Do not update certain indices
            if B[index] and index not in indices_to_keep:

                # Test to see if index is a simple point
                update, n_in = simple_test(faces, index, B, 0, neighbor_lists)

                # If a simple point, remove and run again
                if update and n_in > 1:
                    B[index] = 0
                    exist_simple = True

    return B
