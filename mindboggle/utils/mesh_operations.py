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

    n_vertices: integer
        number of vertices on the mesh

    Returns
    -------
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> n_vertices = 5
    >>> find_neighbors(faces, n_vertices)
        [[1, 2, 3, 4], [0, 2, 4, 3], [0, 1, 3], [0, 2, 4, 1], [0, 3, 1]]

    >>> # Real example:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> points, faces, scalars, n_vertices = load_scalars(depth_file, False)
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> # Write results to vtk file and view with mayavi2:
    >>> index = 0
    >>> IDs = -1 * np.ones(len(points))
    >>> IDs[index] = 1
    >>> IDs[neighbor_lists[index]] = 2
    >>> rewrite_scalar_lists(depth_file, 'test_find_neighbors.vtk',
    >>>                      [IDs], ['neighbors'], IDs)
    >>> os.system('mayavi2 -m Surface -d test_find_neighbors.vtk &')

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
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Examples
    --------
    >>> from mindboggle.utils.mesh_operations import find_neighbors_vertex
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

#-----------------------------------------------------------------------------
# find all triangle faces centered at each node on the mesh
#-----------------------------------------------------------------------------
def find_faces_at_vertices(faces, n_vertices):
    """
    For each vertex, find all faces containing this vertex.
    Note: faces do not have to be triangles.

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero

    n_vertices: integer
        number of vertices on the mesh

    Returns
    --------
    faces_at_vertex : list of lists of integers
        faces_at_vertices[i] is a list of faces that contain the i-th vertex

    Examples
    --------
    >>> # Simple example:
    >>> from mindboggle.utils.mesh_operations import find_faces_at_vertices
    >>> faces = [[0,1,2],[0,2,3],[0,3,4],[0,1,4],[4,3,1]]
    >>> n_vertices = 5
    >>> find_faces_at_vertices(faces, n_vertices)
        [[0, 1, 2, 3], [0, 3, 4], [0, 1], [1, 2, 4], [2, 3, 4]]

    """
    faces_at_vertices = [[] for i in xrange(n_vertices)]
    for face_id, face in enumerate(faces):
        for vertex in face:
           faces_at_vertices[vertex].append(face_id)

    return faces_at_vertices

#------------------------------------------------------------------------------
# Find special "anchor" points for constructing fundus curves
#------------------------------------------------------------------------------
def find_anchors(points, L, min_directions, min_distance, thr):
    """
    Find special "anchor" points for constructing fundus curves.

    Assign maximum likelihood points as "anchor points"
    while ensuring that the anchor points are not close to one another.

    Parameters
    ----------
    points : numpy array of floats
        coordinates for all vertices
    L : numpy array of integers
        fundus likelihood values for all vertices (default -1)
    min_directions : numpy array of floats
        minimum directions for all vertices
    min_distance : minimum distance
    thr : likelihood threshold

    Returns
    -------
    anchors : list of subset of surface mesh vertex indices

    Examples
    --------
    >>> # Use depth values instead of likelihood values for the example
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> from mindboggle.utils.mesh_operations import find_anchors
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.min.dir.txt')
    >>> points, faces, values, n_vertices = load_scalars(depth_file, True)
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> min_distance = 5
    >>> thr = 0.5
    >>> anchors = find_anchors(points, values, min_directions, min_distance, thr)
    >>> # Write results to vtk file and view with mayavi2:
    >>> IDs = -1 * np.ones(len(min_directions))
    >>> IDs[anchors] = 1
    >>> rewrite_scalar_lists(depth_file, 'test_find_anchors.vtk',
    >>>                      [IDs], ['anchors'], IDs)
    >>> os.system('mayavi2 -m Surface -d test_find_anchors.vtk &')

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
# Propagate labels to segment surface into contiguous regions
#------------------------------------------------------------------------------
def propagate(points, faces, region, seeds, labels,
              max_iters=500, tol=0.001, sigma=10):
    """
    Propagate labels to segment surface into contiguous regions,
    starting from seed vertices.

    Parameters
    ----------
    points : list of lists of three integers
        coordinates for all vertices
    faces : list of lists of three integers
        indices to three vertices per face (indices start from zero)
    region : numpy array of integers
        values > -1 indicate inclusion in a region for all vertices
    seeds : numpy array of integers
        seed numbers for all vertices (default -1 for not a seed)
    labels : numpy array of integers
        label numbers for all vertices, with -1s for unlabeled vertices
    max_iters : integer
        maximum number of iterations to run graph-based learning algorithm
    tol: float
        threshold to assess convergence of the algorithm
    sigma: float
        gaussian kernel parameter

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> import mindboggle.label.rebound as rb
    >>> from mindboggle.utils.mesh_operations import find_neighbors,\
    >>>     inside_faces, propagate, detect_boundaries
    >>> from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> import mindboggle.utils.kernels as kernels
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(data_path, 'results', 'features',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'folds.vtk')
    >>> points, faces, folds, n_vertices = load_scalars(folds_file, True)
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, labels, n_vertices = load_scalars(label_file, True)
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> indices_boundaries, label_pairs, foo = detect_boundaries(range(len(points)),
    >>>     labels, neighbor_lists)
    >>> label_pair_lists = sulcus_boundaries()
    >>> fold_ID = 2
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>> fold_array = -1 * np.ones(len(points))
    >>> fold_array[indices_fold] = 1
    >>> indices_boundaries, label_pairs, foo = detect_boundaries(indices_fold,
    >>>     labels, neighbor_lists)
    >>> seeds = -1 * np.ones(len(points))
    >>> for ilist,label_pair_list in enumerate(label_pair_lists):
    >>>     I = [x for i,x in enumerate(indices_boundaries)
    >>>          if np.sort(label_pairs[i]).tolist() in label_pair_list]
    >>>     seeds[I] = ilist
    >>>
    >>> segments = propagate(points, faces, fold_array, seeds, labels)
    >>>
    >>> # Write results to vtk file and view with mayavi2:
    >>> #rewrite_scalar_lists(label_file, 'test_propagate.vtk',
    >>> #                     [segments.tolist()], ['segments'], segments.tolist())
    >>> indices = [i for i,x in enumerate(segments) if x > -1]
    >>> write_scalar_lists('test_propagate.vtk', points, indices,
    >>>     inside_faces(faces, indices), [segments.tolist()], ['segments'])
    >>> os.system('mayavi2 -m Surface -d test_propagate.vtk &')

    """
    import numpy as np
    from mindboggle.utils.mesh_operations import inside_faces
    import mindboggle.utils.kernels as kernels
    import mindboggle.label.rebound as rb

    indices_region = [i for i,x in enumerate(region) if x > -1]
    local_indices_region = -1 * np.ones(len(region))
    local_indices_region[indices_region] = range(len(indices_region))
    points = np.asarray(points)

    n_sets = len(np.unique([x for x in seeds if x > -1]))
    if n_sets == 1:
        print('Segment {0} vertices from 1 set of seed vertices'.
              format(len(indices_region)))
    else:
        print('Segment {0} vertices from {1} sets of seed vertices'.
              format(len(indices_region), n_sets))

    # Set up rebound Bounds class instance
    B = rb.Bounds()
    B.Indices = local_indices_region
    B.Points = points[indices_region, :]
    B.Faces = inside_faces(faces, indices_region)
    B.Labels = labels[indices_region]
    B.seed_labels = seeds[indices_region]
    B.num_points = len(B.Points)

    # Propagate seed IDs from seeds
    B.graph_based_learning(method='propagate_labels', realign=False,
                           kernel=kernels.rbf_kernel, sigma=sigma,
                           max_iters=max_iters, tol=tol, vis=False)

    # Assign maximum probability seed IDs to each point of region
    max_prob_labels = B.assign_max_prob_label()

    # Return segment IDs in original vertex array
    segments = -1 * np.ones(len(points))
    segments[indices_region] = max_prob_labels

    return segments

#------------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions by seed growing
#------------------------------------------------------------------------------
def segment(vertices_to_segment, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False,
            spread_within_labels=False, labels=[], label_lists=[]):
    """
    Segment vertices of surface into contiguous regions by seed growing,
    starting from zero or more lists of seed vertices.

    Parameters
    ----------
    vertices_to_segment : list or array of integers
        indices to mesh vertices to be segmented
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    min_region_size : integer
        minimum size of segmented set of vertices
    seed_lists : list of lists, or empty list
        each list contains indices to seed vertices to segment vertices_to_segment
    keep_seeding : Boolean
        grow from new seeds even after all seed lists have fully grown
    spread_within_labels : Boolean
        grow seeds only by vertices with labels in the seed labels?
    labels : numpy array of integers (required only if spread_within_labels)
        label numbers for all vertices, with -1s for unlabeled vertices
    label_lists : list of lists of integers (required only if spread_within_labels)
        List of unique labels for each seed list to grow into
        (If empty, set to unique labels for each seed list)

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> # Setup
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh_operations import find_neighbors, segment, detect_boundaries
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> vertices_to_segment = np.where(depths > 0.50)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>>
    >>> # Example 1: with seed lists
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> label_pair_lists = sulcus_boundaries()
    >>> label_lists = [np.unique(np.ravel(x)) for x in label_pair_lists]
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, labels, n_vertices = load_scalars(label_file, True)
    >>> indices_boundaries, label_pairs, foo = detect_boundaries(vertices_to_segment,
    >>>     labels, neighbor_lists)
    >>> seed_lists = []
    >>> for label_pair_list in label_pair_lists:
    >>>     seed_lists.append([x for i,x in enumerate(indices_boundaries)
    >>>         if np.sort(label_pairs[i]).tolist() in label_pair_list])
    >>> #seed_lists = [vertices_to_segment[range(2000)],
    >>> #              vertices_to_segment[range(2000,4000)],
    >>> #              vertices_to_segment[range(10000,12000)]]
    >>>
    >>> sulci = segment(vertices_to_segment, neighbor_lists, 50,
    >>>                 seed_lists, True, True, labels, label_lists)
    >>>
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalar_lists(depth_file, 'test_segment.vtk',
    >>>                      [sulci.tolist()], ['sulci'], sulci)
    >>> os.system('mayavi2 -m Surface -d test_segment.vtk &')
    >>>
    >>> # Example 2: without seed lists
    >>> folds = segment(vertices_to_segment, neighbor_lists)
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalar_lists(depth_file, 'test_segment2.vtk',
    >>>                      [folds], ['folds'], folds)
    >>> os.system('mayavi2 -m Surface -d test_segment2.vtk &')

    """
    import numpy as np

    # If seed_lists is empty, select first vertex from vertices_to_segment
    # (single vertex selection does not affect result -- see below*)
    if len(seed_lists):
        select_single_seed = False
        if len(seed_lists) == 1:
            print('    Segment {0} vertices using seed vertices'.
                  format(len(vertices_to_segment)))
        else:
            print('    Segment {0} vertices from {1} sets of seed vertices'.
                  format(len(vertices_to_segment), len(seed_lists)))
    else:
        select_single_seed = True
        seed_lists = [[vertices_to_segment[0]]]
        print('    Segment {0} vertices, selecting the first vertex as initial seed'.
              format(len(vertices_to_segment)))

    # Initialize variables, including the list of vertex indices for each region,
    # vertex indices for all regions, and Boolean list indicating which regions
    # are not fully grown, max. region ID, number of segments, etc.
    segments = -1 * np.ones(len(neighbor_lists))
    region_lists = [[] for x in seed_lists]
    all_regions = []
    fully_grown = [False for x in seed_lists]
    new_segment_index = 0
    counter = 0

    # If label_lists empty, set to unique labels for each seed list
    if spread_within_labels:
        if not len(label_lists):
            label_lists = []
            for seed_list in seed_lists:
                seed_labels = np.unique([labels[x] for x in seed_list])
                label_lists.append(seed_labels)

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
                if len(vertices_to_segment):

                    # Identify neighbors of seeds
                    neighbors = []
                    [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                    # Select neighbors that have not been previously selected
                    # and are among the vertices to segment
                    seed_list = [x for x in list(set(neighbors))
                                 if x not in all_regions
                                 if x in vertices_to_segment]
                else:
                    seed_list = []

                # If there are seeds remaining
                if len(seed_list):

                    # Select neighbors with the same labels
                    # as the initial seed labels
                    if spread_within_labels:
                        seed_list = [x for x in seed_list
                                     if labels[x] in label_lists[ilist]]

                    # Continue growing seed list
                    seed_lists[ilist] = seed_list

                # If there are no seeds remaining
                else:
                    # Stop growing seed list (see exception below)
                    fully_grown[ilist] = True

                    # If the region size is large enough
                    size_region = len(region_lists[ilist])
                    if size_region >= min_region_size:

                        # Assign ID to segmented region and increment ID
                        if select_single_seed:
                            new_segment_index = counter
                            counter += 1
                        else:
                            new_segment_index = ilist
                        segments[region_lists[ilist]] = new_segment_index

                        # Display current number and size of region
                        if size_region > 1:
                            if len(seed_lists) == 1:
                                if len(vertices_to_segment):
                                    print("      {0} vertices remain".
                                          format(len(vertices_to_segment)))
                            else:
                                print("      Region {0}: {1} vertices ({2} remain)".
                                      format(int(new_segment_index), size_region,
                                             len(vertices_to_segment)))

                    # If selecting a single seed, continue growing
                    # if there are more vertices to segment
                    if select_single_seed:
                        if len(vertices_to_segment) >= min_region_size:
                            fully_grown[0] = False
                            seed_lists[0] = [vertices_to_segment[0]]
                            region_lists[0] = []

    # Keep growing from new seeds even after all seed lists have fully grown
    if keep_seeding and len(vertices_to_segment) >= min_region_size:
        print('    Keep seeding to segment {0} remaining vertices'.
              format(len(vertices_to_segment)))

        # Select first unsegmented vertex as new seed
        seed_list = [vertices_to_segment[0]]

        # Loop until the seed list has grown to its full extent
        new_segment_index = ilist + 1
        region = []
        while len(vertices_to_segment) > min_region_size:

            # Add seeds to region
            region.extend(seed_list)
            all_regions.extend(seed_list)

            # Remove seeds from vertices to segment
            vertices_to_segment = list(frozenset(vertices_to_segment).
                                       difference(seed_list))
            if len(vertices_to_segment):

                # Identify neighbors of seeds
                neighbors = []
                [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                # Select neighbors that have not been previously selected
                # and are among the vertices to segment
                seed_list = [x for x in list(set(neighbors))
                             if x not in all_regions
                             if x in vertices_to_segment]
            else:
                seed_list = []

            print(len(vertices_to_segment), len(seed_list))

            # If there are no seeds remaining
            if not len(seed_list):

                # If the region size is large enough
                size_region = len(region)
                if size_region >= min_region_size:

                    # Assign ID to segmented region and increment ID
                    segments[region] = new_segment_index
                    new_segment_index += 1

                    # Display current number and size of region
                    if size_region > 1:
                        print("      {0} vertices remain".
                              format(len(vertices_to_segment)))

                # Select first unsegmented vertex as new seed
                if len(vertices_to_segment) >= min_region_size:
                    seed_list = [vertices_to_segment[0]]
                    region = []

    return segments

#------------------------------------------------------------------------------
# Fill holes
#------------------------------------------------------------------------------
def label_holes(holes, hole_numbers, regions, neighbor_lists):
    """
    Fill holes in regions on a surface mesh.

    Parameters
    ----------
    holes : numpy array of integers
        hole numbers for all vertices (default -1)
    hole_numbers : list or numpy array of integers
        unique hole numbers (no -1)
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

    # Make sure arguments are numpy arrays
    if type(regions) != np.ndarray:
        regions = np.array(regions)

    # Identify the vertices for each hole
    for n_hole in hole_numbers:
        I = np.where(holes == n_hole)[0]

        # Identify neighbors to these vertices
        N=[]; [N.extend(neighbor_lists[i]) for i in I]
        if len(N):

            # Assign the hole the maximum region ID number of its neighbors
            regions[I] = max([regions[x] for x in N])

    return regions

def fill_holes(regions, neighbor_lists):
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

    Returns
    -------
    regions : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh_operations import find_neighbors, segment
    >>> from mindboggle.utils.mesh_operations import inside_faces, fill_holes
    >>> from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> # Select one sulcus
    >>> sulci_file = os.path.join(data_path, 'results', 'features',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    >>> points, faces, sulci, n_vertices = load_scalars(sulci_file, True)
    >>> #neighbor_lists = find_neighbors(faces, n_vertices)
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> n_sulcus = 0
    >>> sulci[sulci != n_sulcus] = -1
    >>> # Make hole in sulcus in case there isn't one
    >>> I = np.where(sulci==n_sulcus)[0]
    >>> index = I[np.round(len(I)/2)]
    >>> N = neighbor_lists[index]
    >>> for n in N:
    >>>     if any(sulci[neighbor_lists[n]] == -1):
    >>>         print("Select a different index")
    >>>         break
    >>> neighbor_lists[index] = []
    >>> for n in N:
    >>>     neighbor_lists[n] = []
    >>> sulci[index] = -1
    >>> sulci[N] = -1
    >>> # Write hole to vtk file and view with mayavi2:
    >>> indices = [i for i,x in enumerate(sulci) if x > -1]
    >>> write_scalar_lists('test_hole.vtk', points, indices,
    >>>     inside_faces(faces, indices), [sulci], ['holes'])
    >>> os.system('mayavi2 -m Surface -d test_hole.vtk &')
    >>> # Fill hole
    >>> regions = fill_holes(sulci, neighbor_lists)
    >>> # Write results to vtk file and view with mayavi2:
    >>> indices = [i for i,x in enumerate(regions) if x > -1]
    >>> write_scalar_lists('test_fill_holes.vtk', points, indices,
    >>>     inside_faces(faces, indices), [regions], ['regions'])
    >>> os.system('mayavi2 -m Surface -d test_fill_holes.vtk &')

    """
    import numpy as np

    #--------------------------------------------------------------------------
    # Find boundaries to holes
    #--------------------------------------------------------------------------
    hole_boundaries = -1 * np.ones(len(regions))

    # Identify vertices for each region
    region_numbers = [x for x in np.unique(regions) if x > -1]
    count = 0
    for n_region in region_numbers:
        region_indices = np.where(regions == n_region)[0]

        # Identify neighbors to these vertices and their neighbors
        N = []
        [N.extend(neighbor_lists[x]) for x in region_indices]
        N = np.unique(N)
        N = [x for x in N if x not in region_indices]
        N2 = []
        [N2.extend(neighbor_lists[x]) for x in N if x not in region_indices]
        N.extend(N2)
        N = np.unique(N)
        if len(N):

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

        # Grow seeds from hole boundaries to fill holes (background value -1)
        seed_lists = []
        for n_hole in hole_numbers:
            I = np.where(hole_boundaries == n_hole)[0]
            seed_lists.append(I)
        background = [i for i,x in enumerate(regions) if x == -1]
        holes = segment(background, neighbor_lists, 1, seed_lists)

        # Label the vertices for each hole by surrounding region number
        regions = label_holes(holes, hole_numbers, regions, neighbor_lists)

    return regions

#------------------------------------------------------------------------------
# Test for simple points
#------------------------------------------------------------------------------
def simple_test(index, values, neighbor_lists):
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

    # Make sure arguments are numpy arrays
    if type(values) != np.ndarray:
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
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    >>> from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    >>> from mindboggle.extract.extract_fundi import compute_likelihood, connect_points
    >>> from mindboggle.utils.mesh_operations import find_anchors, skeletonize, extract_endpoints
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.avg.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.min.dir.txt')
    >>> sulci_file = os.path.join(data_path, 'results', 'features',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    >>> points, faces, sulci, n_vertices = load_scalars(sulci_file, True)
    >>> points, faces, mean_curvatures, n_vertices = load_scalars(mean_curvature_file,
    >>>                                                          True)
    >>> points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> n_sulcus = max(sulci)
    >>> sulci[sulci != n_sulcus] = -1
    >>> indices_sulcus = [i for i,x in enumerate(sulci) if x == n_sulcus]
    >>> sulcus_likelihoods = compute_likelihood(depths[indices_sulcus],
    >>>                                       mean_curvatures[indices_sulcus])
    >>> thr = 0.5
    >>> L = np.zeros(len(points))
    >>> L[indices_sulcus] = sulcus_likelihoods
    >>> sulcus_indices_anchors = find_anchors(points[indices_sulcus, :],
    >>>     sulcus_likelihoods, min_directions[indices_sulcus], 5, thr)
    >>> indices_anchors = [indices_sulcus[x] for x in sulcus_indices_anchors]
    >>> L[L > thr] = 1
    >>> L[L <= thr] = 0
    >>> anchor_skeleton = skeletonize(L, indices_anchors, neighbor_lists)
    >>> indices_skeleton = [i for i,x in enumerate(anchor_skeleton) if x > 0]
    >>> indices_endpoints = extract_endpoints(indices_skeleton, neighbor_lists)
    >>> indices_endpoints = [x for x in indices_endpoints if x in indices_anchors]
    >>> skeleton = skeletonize(L, indices_endpoints, neighbor_lists)
    >>> indices_endpoint_skel = [x for x in skeleton if x > 0]
    >>> # Write results to vtk file and view with mayavi2:
    >>> skeleton[indices_anchors] = 3
    >>> skeleton[indices_endpoints] = 5
    >>> indices = [i for i,x in enumerate(sulci) if x > -1]
    >>> write_scalar_lists('test_skeletonize.vtk', points, indices,
    >>>     inside_faces(faces, indices), [skeleton], ['skeleton'])
    >>> os.system('mayavi2 -m Surface -d test_skeletonize.vtk &')

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
                update, n_in = simple_test(index, binary_array, neighbor_lists)

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
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.utils.mesh_operations import detect_boundaries, extract_endpoints
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> label_pair_lists = sulcus_boundaries()
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, labels, n_vertices = load_scalars(label_file, True)
    >>> label_indices = [i for i,x in enumerate(labels)
    >>>                  if x in label_pair_lists[0][0]]
    >>> indices_boundaries, label_pairs, foo = detect_boundaries(label_indices,
    >>>                                               labels, neighbor_lists)
    >>> indices_endpoints = extract_endpoints(indices_boundaries, neighbor_lists)
    >>> # Write results to vtk file and view with mayavi2:
    >>> end_IDs = -1 * np.ones(len(points))
    >>> end_IDs[indices_boundaries] = 1
    >>> end_IDs[indices_endpoints] = 2
    >>> rewrite_scalar_lists(depth_file, 'test_extract_endpoints.vtk',
    >>>                      [end_IDs], ['endpoints'], end_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_endpoints.vtk &')

    """

    # Find vertices with only one connected neighbor
    indices_endpoints = []
    for index in indices_skeleton:
        if len([x for x in neighbor_lists[index] if x in indices_skeleton]) == 1:
            indices_endpoints.append(index)

    return indices_endpoints

#------------------------------------------------------------------------------
# Filter faces
#------------------------------------------------------------------------------
def inside_faces(faces, indices):
    """
    Remove surface mesh faces whose three vertices are not all in "indices"

    Parameters
    ----------
    faces : list of lists of three integers
        the integers for each face are indices to vertices, starting from zero
    indices : vertex indices to mesh

    Returns
    -------
    faces : reduced array of faces

    Examples
    --------
    >>> from mindboggle.utils.mesh_operations import inside_faces
    >>> faces = [[1,2,3], [2,3,7], [4,7,8], [3,2,5]]
    >>> indices = [0,1,2,3,4,5]
    >>> inside_faces(faces, indices)
      Reduced 4 to 2 triangular faces.
      array([[1, 2, 3],
             [3, 2, 5]])

    """
    import numpy as np

    len_faces = len(faces)
    fs = frozenset(indices)
    faces = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces = np.reshape(np.ravel(faces), (-1, 3))
    if len(faces) < len_faces:
        print('Reduced {0} to {1} triangular faces'.format(len_faces, len(faces)))

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
        label numbers for all vertices, with -1s for unlabeled vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    boundary_indices : list of integers
        indices to label boundary vertices
    boundary_label_pairs : list of lists of sorted pairs of integers
        label pairs
    unique_boundary_label_pairs : list of pairs of integers
        unique label pairs

    Examples
    --------
    >>> # Small example:
    >>> from mindboggle.utils.mesh_operations import detect_boundaries
    >>> neighbor_lists = [[1,2,3], [0,0,8,0,8], [2], [4,7,4], [3,2,3]]
    >>> labels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    >>> region = [0,1,2,4,5,8,9]
    >>> detect_boundaries(region, labels, neighbor_lists)
      ([1, 4], [[10, 90], [40, 30]])

    >>> Real example -- extract sulcus label boundaries:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.utils.mesh_operations import detect_boundaries
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, labels, n_vertices = load_scalars(label_file, True)
    >>> neighbor_lists = find_neighbors(faces, n_vertices)
    >>> indices_boundaries, label_pairs, foo = detect_boundaries(range(len(points)),
    >>>     labels, neighbor_lists)
    >>> # Write results to vtk file and view with mayavi2:
    >>> IDs = -1 * np.ones(len(points))
    >>> IDs[indices_boundaries] = 1
    >>> rewrite_scalar_lists(label_file, 'test_detect_boundaries.vtk',
    >>>                      [IDs], ['boundaries'], IDs)
    >>> os.system('mayavi2 -m Surface -d test_detect_boundaries.vtk &')

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

#------------------------------------------------------------------------------
# Compute distance
#------------------------------------------------------------------------------
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

    Examples
    --------
    >>> from mindboggle.utils.mesh_operations import compute_distance
    >>> point = [1,2,3]
    >>> points = [[10,2.0,3], [0,1.5,2]]
    >>> compute_distance(point, points)
      (1.5, 1)

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

# propagate() example
if __name__ == "__main__" :
    import os
    import numpy as np
    import mindboggle.label.rebound as rb
    from mindboggle.utils.mesh_operations import find_neighbors, inside_faces, \
        detect_boundaries
    from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    import mindboggle.utils.kernels as kernels

    data_path = os.environ['MINDBOGGLE_DATA']
    label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
                 'label', 'lh.labels.DKT25.manual.vtk')
    points, faces, labels, n_vertices = load_scalars(label_file, True)

    neighbor_lists = find_neighbors(faces, n_vertices)

    indices_boundaries, label_pairs, foo = detect_boundaries(range(len(points)),
        labels, neighbor_lists)

    sulcuss_file = os.path.join(data_path, 'results', 'features',
                 '_hemi_lh_subject_MMRR-21-1', 'sulcuss.vtk')
    points, faces, folds, n_vertices = load_scalars(folds_file, True)

    label_pair_lists = sulcus_boundaries()

    fold_ID = 1
    indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    fold_array = np.zeros(len(points))
    fold_array[indices_fold] = 1
    indices_boundaries, label_pairs, foo = detect_boundaries(indices_fold,
        labels, neighbor_lists)

    seeds = -1 * np.ones(len(points))
    for ilist,label_pair_list in enumerate(label_pair_lists):
        I = [x for i,x in enumerate(indices_boundaries)
             if np.sort(label_pairs[i]).tolist() in label_pair_list]
        seeds[I] = ilist

    segments = propagate(points, faces, fold_array, seeds, labels)

    # Write results to vtk file and view with mayavi2:
    rewrite_scalar_lists(label_file, 'test_propagate.vtk',
                         [segments], ['segments'], segments)
    os.system('mayavi2 -m Surface -d test_propagate.vtk &')
