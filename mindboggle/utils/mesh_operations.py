#!/usr/bin/python
"""
Operations on surface mesh vertices.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Bao  (forrest.bao@gmail.com)

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import numpy as np
#from operator import itemgetter

#------------------------------------------------------------------------------
# Find all neighbors from faces
#------------------------------------------------------------------------------
def find_neighbors(faces, n_vertices):
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
    >>> from utils.mesh_operations import find_neighbors
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> n_vertices = 5
    >>> find_neighbors(faces, n_vertices)
        [[1, 2, 3, 4], [0, 2, 4, 3], [1, 0, 3], [2, 0, 4, 1], [0, 3, 1]]

    >>> from utils.mesh_operations import find_neighbors
    >>> from utils.io_vtk import load_scalar
    >>> points, faces, scalars, n_vertices = load_scalar('lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors(faces, n_vertices)

    """

    neighbor_lists = [[] for x in xrange(n_vertices)]

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
    neighbor_list : list of integers
        indices of neighboring vertices

    Example
    -------
    >>> from utils.mesh_operations import find_neighbors_vertex
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4]]
    >>> index = 1
    >>> find_neighbors_vertex(faces, index)
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

#------------------------------------------------------------------------------
# Find anchor points
#------------------------------------------------------------------------------
def find_anchors(points, L, min_directions, min_distance, thr):
    """
    Find anchor points.

    Assign maximum likelihood points as "anchor points"
    while ensuring that the anchor points are not close to one another.
    Anchor points are used to construct curves.

    Note: only the sulcus end points are strictly necessary
          (so that the fundus doesn't shrink)

    Parameters
    ----------
    points : [#points x 3] numpy array
    L : fundus likelihood values: [#points x 1] numpy array
    min_directions : [#points x 1] numpy array
    min_distance : minimum distance
    thr : likelihood threshold

    Returns
    -------
    anchors : list of subset of surface mesh vertex indices

    Example
    -------
    >>> import numpy as np
    >>> from utils.io_vtk import load_scalar
    >>> from utils.mesh_operations import find_anchors
    >>> min_curvature_vector_file = '/desk/output/results/measures/' + \
    >>>     '_hemi_lh_subject_MMRR-21-1/lh.pial.curv.min.dir.txt'
    >>> values_file = '/desk/output/results/measures/' + \
    >>>     '_hemi_rh_subject_MMRR-21-1/lh.pial.depth.vtk'
    >>> points, faces, values, n_vertices = load_scalar(values_file, 1)
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> min_distance = 5
    >>> thr = 0.5
    >>> anchors = find_anchors(range(len(values)),
    >>>                        values, min_directions, min_distance, thr)
    >>> # Write results to vtk file:
    >>> from utils.io_vtk import rewrite_scalars
    >>> LUT = np.zeros(len(min_directions))
    >>> LUT[anchors] = 1
    >>> rewrite_scalars(values_file, 'test_find_anchors.vtk', LUT)

    """
    import numpy as np
    from operator import itemgetter

    # Make sure arguments are numpy arrays
    if type(points) != np.ndarray:
        points = np.array(points)
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
            D = np.linalg.norm(points[anchors[i], :] - points[imax, :])

            # If distance less than threshold, consider the point found
            if D < min_distance:
                found = 1
            # Compute directional distance between points if they are close
            elif D < max_distance:
                dirV = np.dot(points[anchors[i], :] - points[imax, :],
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

#------------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions
#------------------------------------------------------------------------------
def segment(vertices_to_segment, neighbor_lists, seed_lists=[], min_region_size=1,
            spread_same_labels=False, labels=[], label_pair_lists=[]):
    """
    Segment vertices of surface into contiguous regions,
    starting from zero or more lists of seed vertices.

    Parameters
    ----------
    vertices_to_segment : list or array of integers
        indices to subset of mesh vertices to be segmented
    neighbor_lists : list of lists of integers (#lists = #vertices in mesh)
        each list contains indices to neighboring vertices
    seed_lists : list of lists, or empty list
        each list contains indices to seed vertices to segment vertices_to_segment
    min_region_size : minimum size of segmented set of vertices
    spread_same_labels : Boolean
        grow seeds only by vertices with labels in the seed labels or not
    labels : numpy array of integers (optional)
        labels, one per vertex in the mesh
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus

    Returns
    -------
    segments : numpy array [total #vertices in mesh x 1]
        segment indices for regions (default -1)
    n_segments : #segments

    Example
    -------
    >>> import numpy as np
    >>> from utils.mesh_operations import find_neighbors, segment
    >>> from utils.io_vtk import load_scalar, rewrite_scalars
    >>> from info.sulcus_boundaries import sulcus_boundaries
    >>> depth_file = '/desk/output/results/measures/_hemi_lh_subject_MMRR-21-1/lh.pial.depth.vtk'
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, True)
    >>> vertices_to_segment = np.where(depths > .25)[0]
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> folds, n_folds = segment(vertices_to_segment, neighbor_lists,
    >>>     seed_lists=[], min_region_size=1,
    >>>     spread_same_labels=False, labels=[], label_pair_lists=[])
    >>> # Write results to vtk file:
    >>> rewrite_scalars(depth_file, 'test_segment1.vtk', folds, folds)

    >>> seed_lists = [vertices_to_segment[range(100)],
                      vertices_to_segment[range(100,200)],
                      vertices_to_segment[range(200,300)]]
    >>> label_file = '/Applications/freesurfer/subjects/MMRR-21-1/label/lh.labels.DKT25.manual.vtk'
    >>> points, faces, labels, n_vertices = load_scalar(label_file, True)
    >>> label_pair_lists = sulcus_boundaries()
    >>> folds, n_folds = segment(vertices_to_segment, neighbor_lists, seed_lists,
    >>>     50, True, labels, label_pair_lists)
    >>> # Write results to vtk file:
    >>> rewrite_scalars(depth_file, 'test_segment2.vtk', folds, folds)

    """
    import numpy as np

    # If seed_lists is empty, select first vertex from vertices_to_segment
    # (single vertex selection does not affect result -- see below*)
    print('    Segment {0} vertices...'.format(len(vertices_to_segment)))
    if len(seed_lists):
        select_single_seed = False
        print('    Select from {0} sets of seed vertices'.format(len(seed_lists)))
    else:
        select_single_seed = True
        seed_lists = [[vertices_to_segment[0]]]
        print('    Select first vertex as initial seed')

    # Initialize variables, including the list of vertex indices for each region,
    # vertex indices for all regions, and Boolean list indicating which regions
    # are not fully grown, max. region ID, number of segments, etc.
    segments = -1 * np.ones(len(neighbor_lists))
    region_lists = [[] for x in seed_lists]
    all_regions = []
    fully_grown = [False for x in seed_lists]
    new_segment_index = 0
    counter = 0

    # Prepare list of all unique sorted label pairs in the labeling protocol
    label_lists = [np.ravel(x) for x in label_pair_lists]

    # Loop until all of the seed lists have grown to their full extent
    while not all(fully_grown):
        # Loop through seed lists over and over again
        for ilist, seed_list in enumerate(seed_lists):
            # If seed list empty
            if not len(seed_list):
                fully_grown[ilist] = True
            # If seed list not fully grown
            if not fully_grown[ilist]:

                # Add seeds to region
                region_lists[ilist].extend(seed_list)
                all_regions.extend(seed_list)

                # Remove seeds from vertices to segment
                vertices_to_segment = list(frozenset(vertices_to_segment).
                                        difference(seed_list))

                # Identify neighbors of seeds
                neighbors = []
                [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                # Select neighbors that have not been previously selected
                # and are among the vertices to segment
                seed_list = [x for x in list(set(neighbors))
                             if x not in all_regions
                             if x in vertices_to_segment]

                # Select neighbors with the same labels
                # as the initial seed label pairs
                if len(seed_list):
                    if spread_same_labels and not select_single_seed:
                        print(labels[list(set(neighbors))])
                        print(label_lists[ilist])
                        seed_list = [x for x in seed_list
                                     if labels[x] in label_lists[ilist]]

                print(ilist, len(seed_list), len(region_lists[ilist]))

                # Continue growing seed list if there are seeds remaining
                if len(seed_list):
                    seed_lists[ilist] = seed_list

                # If there are no seeds remaining
                else:
                    # Stop growing seed list (see exception below)
                    fully_grown[ilist] = True

                    # If the region size is large enough
                    size_region = len(region_lists[ilist])
                    if size_region >= min_region_size:

                        # Assign ID to segmented region
                        if select_single_seed:
                            new_segment_index = counter
                            counter += 1
                        else:
                            new_segment_index = ilist
                        segments[region_lists[ilist]] = new_segment_index

                        # Display current number and size of region
                        if size_region > 1:
                            print("    Region {0}: {1} vertices "
                                  "({2} remaining)".
                                  format(new_segment_index, size_region,
                                         len(vertices_to_segment)))

                    # If selecting a single seed, continue growing
                    # if there are more vertices to segment
                    if select_single_seed:
                        if len(vertices_to_segment) >= min_region_size:
                            fully_grown[0] = False
                            seed_lists[0] = [vertices_to_segment[0]]
                            region_lists[0] = []

    n_segments = new_segment_index + 1
    return segments, n_segments

#------------------------------------------------------------------------------
# Fill holes
#------------------------------------------------------------------------------
def fill_holes(regions, holes, n_holes, neighbor_lists):
    """
    Fill holes in surface mesh regions.

    Parameters
    ----------
    regions : numpy array [total #vertices in mesh x 1]
        indices for regions (default -1)
    holes : [#vertices x 1] numpy array
    n_holes : [#vertices x 1] numpy array
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices

    Returns
    -------
    regions : [#vertices x 1] numpy array of integers
        region ID numbers

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if type(regions) != np.ndarray:
        regions = np.array(regions)

    # Identify the vertices for each hole
    for n_hole in range(n_holes):
        I = np.where(holes == n_hole)[0]

        # Identify neighbors to these vertices
        N=[]; [N.extend(neighbor_lists[i]) for i in I]
        if len(N):

            # Assign the hole the maximum region ID number of its neighbors
            regions[I] = max([regions[x] for x in N])

    return regions

#------------------------------------------------------------------------------
# Test for simple points
#------------------------------------------------------------------------------
def simple_test(index, values, thr, neighbor_lists):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    Parameters
    ----------
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
                        if len(frozenset(N[i]).intersection(N[j])) > 0:
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

    Parameters
    ----------
    binary_array : binary [#vertices x 1] numpy array
    indices_to_keep : indices to retain
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    binary_array : binary skeleton: numpy array

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if type(binary_array) != np.ndarray:
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
                update, n_in = simple_test(index, binary_array, 0, neighbor_lists)

                # If a simple point, remove and run again
                if update and n_in > 1:
                    binary_array[index] = 0
                    exist_simple = True

    return binary_array


def inside_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices"

    Parameters
    ----------
    faces : triangular surface mesh vertex indices [#faces x 3]
    indices : vertex indices to mesh

    Returns
    -------
    faces : reduced array of faces

    """
    import numpy as np

    len_faces = len(faces)
    fs = frozenset(indices)
    faces = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces = np.reshape(np.ravel(faces), (-1, 3))
    print('  Reduced {0} to {1} triangular faces.'.format(len_faces, len(faces)))

    return faces

#------------------------------------------------------------------------------
# Detect label boundaries
#------------------------------------------------------------------------------
def detect_boundaries(region, labels, neighbor_lists):
    """
    Detect the label boundaries in a collection of vertices such as a region.

    Label boundaries are the set of all vertices
    whose neighbors do not share the same label.

    Parameters
    ----------
    region : list of integers
        indices to vertices in a region (any given collection of vertices)
    labels : numpy array of integers
        labels for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices

    Returns
    -------
    boundary_indices : list of integers
        indices to label boundary vertices
    boundary_label_pairs : list of lists of sorted pairs of integers
        label pairs

    Example
    -------
    >>> from utils.mesh_operations import detect_boundaries
    >>> neighbor_lists = [[1,2,3], [0,0,8,0,8], [2], [4,7,4], [3,2,3]]
    >>> labels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    >>> region = [0,1,2,4,5,8,9]
    >>> detect_boundaries(region, labels, neighbor_lists)
      ([1, 4], [[10, 90], [40, 30]])

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if type(labels) != np.ndarray:
        labels = np.array(labels)

    # Construct a list of labels corresponding to the neighbor lists
    label_lists = [list(set(labels[lst])) for lst in neighbor_lists]

    # Find indices to sets of two labels
    boundary_indices = [i for i,x in enumerate(label_lists)
                        if len(set(x)) == 2
                        if i in region]
    boundary_label_pairs = [x for i,x in enumerate(label_lists)
                            if len(set(x)) == 2
                            if i in region]

    unique_boundary_label_pairs = []
    for pair in boundary_label_pairs:
        if np.sort(pair).tolist() not in unique_boundary_label_pairs:
            unique_boundary_label_pairs.append(np.sort(pair).tolist())

    return boundary_indices, boundary_label_pairs, unique_boundary_label_pairs


def compute_distance(point, points):
    """
    Estimate the Euclidean distance from one point to a second (set) of points.

    Parameters
    ----------
    point : list of three floats
        coordinates for a single point
    points : list with one or more lists of three floats
        coordinates for a second point (or multiple points)

    Returns
    -------
    min_distance : float
        Euclidean distance between two points,
        or the minimum distance between a point and a set of points
    min_index : int
        index of closest of the points (zero if only one)

    """
    import numpy as np

    # If points is a single point
    if np.ndim(points) == 1:
        return np.sqrt((point[0] - points[0]) ** 2 + \
                       (point[1] - points[1]) ** 2 + \
                       (point[2] - points[2]) ** 2), 0

    # If points is a set of multiple points
    elif np.ndim(points) == 2:
        min_distance = np.Inf
        min_index = 0
        for index, point2 in enumerate(points):
            distance = np.sqrt((point[0] - point2[0]) ** 2 + \
                               (point[1] - point2[1]) ** 2 + \
                               (point[2] - point2[2]) ** 2)
            if distance < min_distance:
                min_distance = distance
                min_index = index
        return min_distance, min_index

    # Else return None
    else:
        return None, None
