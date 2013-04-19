#!/usr/bin/env python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame, 2012-2013  .  yrjo.hame@gmail.com
Arno Klein, 2012-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#==================
# Extract all fundi
#==================
def extract_fundi(folds_or_file, depth_file, likelihoods_or_file, save_file=False):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Steps ::

        1. Find fundus endpoints from likelihood and
           minimum distance values and minimum directions.

        2. Connect fundus endpoints and extract fundi.

    Parameters
    ----------
    folds_or_file : list or string
        fold number for each vertex or name of VTK file containing folds scalars
    depth_file :  string
        surface mesh file in VTK format with scalar rescaled depth values
    likelihoods_or_file : list or string
        fundus likelihood values or name of VTK file with the scalar values
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundi : list of integers
        fundus numbers for all vertices (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundus_endpoints : list of integers
        indices to fundus anchor vertices
    fundi_file : string (if save_file)
        name of output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')
    >>> likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> single_fold = False
    >>> if single_fold:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> else:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> #
    >>> fundi, n_fundi, endpoints, fundi_file = extract_fundi(folds, depth_file,
    >>>     likelihoods_or_file, save_file=True)
    >>> #
    >>> # View:
    >>> plot_vtk('fundi.vtk')

    """
    import os
    import sys
    import numpy as np
    from time import time

    from mindboggle.features.fundi import find_endpoints, connect_points
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Load fold numbers:
    if isinstance(folds_or_file, str):
        folds, name = read_scalars(folds_or_file)
    elif isinstance(folds_or_file, list):
        folds = folds_or_file
    elif isinstance(folds_or_file, np.ndarray):
        folds = folds_or_file.tolist()
    else:
        sys.error('folds_or_file is not a string, list, or array.')

    # Load likelihood values:
    if isinstance(likelihoods_or_file, str):
        likelihoods, name = read_scalars(likelihoods_or_file, True, True)
    elif isinstance(likelihoods_or_file, list):
        likelihoods = np.array(likelihoods_or_file)
    elif isinstance(likelihoods_or_file, np.ndarray):
        likelihoods = likelihoods_or_file

    # Initialize variables:
    t1 = time()
    count = 0
    neighbor_lists = find_neighbors_from_file(depth_file)
    depths, name = read_scalars(depth_file)
    npoints = len(folds)
    Z = np.zeros(npoints)
    fundi = -1 * np.ones(npoints)
    fundus_endpoints = []

    # For each fold region...
    unique_fold_IDs = [x for x in np.unique(folds) if x > -1]
    n_folds = len(unique_fold_IDs)
    print("Extract a fundus from each of {0} regions...".format(n_folds))
    for fold_ID in unique_fold_IDs[2:5]:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Region {0}:'.format(int(fold_ID)))

            # Find fundus points
            indices_endpoints, tracks = find_endpoints(indices_fold,
                neighbor_lists, likelihoods, depths, min_edges=5,
                use_threshold=True)
            n_endpoints = len(indices_endpoints)
            if n_endpoints > 1:
                fundus_endpoints.extend(indices_endpoints)
                t2 = time()

                # Connect fundus points and extract fundi:
                print('    Connect {0} fundus points...'.format(n_endpoints))
                B = connect_points(indices_endpoints, indices_fold,
                                   likelihoods, neighbor_lists)
                indices_skeleton = [i for i,x in enumerate(B) if x != -1]
                if len(indices_skeleton) > 1:
                    fundi[indices_skeleton] = fold_ID
                    count += 1
                    print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                          format(n_endpoints, time() - t2))

    n_fundi = count
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.format(n_fundi, time() - t1))

    #---------------------------------------------------------------------------
    # Return fundi, number of fundi, fundus points, and file name:
    #---------------------------------------------------------------------------
    fundi = fundi.tolist()

    if save_file:
        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(depth_file, fundi_file, fundi, 'fundi', folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundus_endpoints, fundi_file

def track_to_border(seed, segments, neighbor_lists, values, borders):
    """
    Build a track from a point through concentric segments.

    This function builds a track from an initial point through concentric
    segments along high-value vertices of a surface mesh.

    Parameters
    ----------
    seed : integer
        index to initial seed vertex from which to grow a track
    segments : list of lists of integers
        indices to vertices for each concentric segment
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    values : numpy array of floats
        values for all vertices that help to guide a track

    Returns
    -------
    track : list of integers
        indices of ordered vertices for a single track

    Examples
    --------
    >>> # Track from deepest point in a fold to its boundary:
    >>> import os
    >>> import pickle
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.labels.label import extract_borders
    >>> from mindboggle.features.fundi import track_to_border
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != -1]
    >>> # Start from the boundary of a thresholded indices:
    >>> use_threshold = True
    >>> if use_threshold:
    >>>     segments_file = os.path.join(path, 'tests', 'segments_fold11.pkl')
    >>>     segments = pickle.load(open(segments_file, 'rb'))
    >>>     seed = segments[0][np.argmax(values[segments[0]])]
    >>> # Or from the maximum value point:
    >>> else:
    >>>     segments_file = os.path.join(path, 'tests', 'segments_likelihood_fold11.pkl')
    >>>     segments = pickle.load(open(segments_file, 'rb'))
    >>>     seed = indices[np.argmax(values[indices])]
    >>> # Extract boundary:
    >>> D = np.ones(len(values))
    >>> D[indices] = 2
    >>> borders, foo1, foo2 = extract_borders(range(len(values)), D, neighbor_lists)
    >>> #
    >>> track = track_to_border(seed, segments, neighbor_lists, values, borders)
    >>> #
    >>> # View:
    >>> T = -1 * np.ones(len(values))
    >>> T[track] = 1
    >>> rewrite_scalars(depth_file, 'track_to_border.vtk', T, 'track', fold)
    >>> plot_vtk('track_to_border.vtk')

    """
    import numpy as np

    track = []
    for isegment, segment in enumerate(segments):

        # Find the seed's neighborhood N in the segment:
        N = neighbor_lists[seed]
        N = [x for x in N if values[x] != -1]
        N_segment = list(frozenset(N).intersection(segment))
        if N:

            # Add the neighborhood vertex with the maximum value to the track:
            if N_segment:
                seed = N_segment[np.argmax(values[N_segment])]
                track.append(seed)

                # If the track has run into the region's border, return the track:
                #if list(frozenset(N_segment).intersection(borders)):
                if seed in borders:
                    return track

            # If there is no neighbor in the new segment,
            # back up to the previous segment:
            elif isegment > 0:
                bridge = []
                max_bridge = 0
                N_previous = list(frozenset(N).intersection(segments[isegment-1]))
                for Np in N_previous:
                    N_next = list(frozenset(neighbor_lists[Np]).intersection(segment))
                    if N_next:
                        if np.max(values[N_next]) > max_bridge:
                            seed = N_next[np.argmax(values[N_next])]
                            bridge = [Np, seed]
                            max_bridge = np.max(values[N_next])
                if bridge:
                    track.extend(bridge)

                    # If the track has run into the region's border, return the track:
                    if seed in borders:
                        return track

        # If there is no neighborhood for the seed, return the track:
        else:
            return None

    # If the track remains empty or does not reach the border, return the track:
    return None

#------------------------------------------------------------------------------
# Find endpoints
#------------------------------------------------------------------------------
def find_endpoints(indices, neighbor_lists, values, values_seeding, min_edges=5,
                   use_threshold=True):
    """
    Find endpoints in a region of connected vertices.

    These points are intended to serve as endpoints of fundus curves
    running along high-likelihood paths within a region (fold).
    This algorithm iteratively propagates paths from a high-likelihood boundary
    within a region of a surface mesh to the boundary of the region.

    Parameters
    ----------
    indices : list of integers
        indices of the vertices to segment (such as a fold in a surface mesh)
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    values : numpy array of floats
        values for all vertices (e.g., fundus likelihood values)
    values_seeding : numpy array of floats
        values for all vertices to threshold for initial seeds
    min_edges : integer
        minimum number of edges between endpoint vertices
    use_threshold : Boolean
        initialize seeds with thresholded vertices?

    Returns
    -------
    indices_endpoints : list of integers
        indices of surface mesh vertices that are endpoints
    tracks : list of lists of integers
        indices to track vertices

    Examples
    --------
    >>> # Setup:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    >>> from mindboggle.features.fundi import find_endpoints
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> values_seeding_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> values_seeding, name = read_scalars(values_seeding_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> min_size = 50
    >>> min_edges = 5
    >>> use_threshold = True
    >>> #
    >>> #-----------------------------------------------------------------------
    >>> # Find endpoints on a single fold:
    >>> #-----------------------------------------------------------------------
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != -1]
    >>> #
    >>> indices_endpoints, tracks = find_endpoints(indices, neighbor_lists, \
    >>>     values, values_seeding, min_edges, use_threshold)
    >>> #
    >>> # Write results to VTK file and view:
    >>> indices_tracks = [x for lst in tracks for x in lst]
    >>> values[indices_tracks] = max(values) + 0.1
    >>> values[indices_endpoints] = max(values) + 0.2
    >>> rewrite_scalars(depth_file, 'find_endpoints.vtk', \
    >>>                 values, 'endpoints_on_values_in_fold', fold)
    >>> plot_vtk('find_endpoints.vtk')
    >>> #-----------------------------------------------------------------------
    >>> # Find endpoints on every fold in a hemisphere:
    >>> #-----------------------------------------------------------------------
    >>> plot_each_fold = 1#False
    >>> folds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> fold_numbers = [x for x in np.unique(folds) if x != -1]
    >>> nfolds = len(fold_numbers)
    >>> endpoints = []
    >>> for ifold, fold_number in enumerate(fold_numbers[0:1]):
    >>>     print('Fold {0} ({1} of {2})'.format(int(fold_number), ifold+1, nfolds))
    >>>     indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>>     if len(indices) > min_size:
    >>>         indices_endpoints, tracks = find_endpoints(indices, neighbor_lists, \
    >>>             values, values_seeding, min_edges, use_threshold)
    >>>         endpoints.extend(indices_endpoints)
    >>>         # Plot each fold:
    >>>         if plot_each_fold:
    >>>             fold = -1 * np.ones(len(values))
    >>>             fold[indices] = 1
    >>>             values[indices_endpoints] = max(values) + 0.1
    >>>             rewrite_scalars(depth_file, 'find_endpoints.vtk',
    >>>                     values, 'endpoints_on_values_in_fold', fold)
    >>>             plot_vtk('find_endpoints.vtk')
    >>> E = -1 * np.ones(len(values))
    >>> E[endpoints] = 1
    >>> #
    >>> # Write results to VTK file and view:
    >>> rewrite_scalars(folds_file, 'find_endpoints.vtk',
    >>>                 E, 'endpoints_on_folds', folds)
    >>> plot_vtk('find_endpoints.vtk')

    """
    import numpy as np

    from mindboggle.labels.label import extract_borders
    from mindboggle.labels.segment import segment, segment_rings
    from mindboggle.features.fundi import track_to_border

    # Initialize R, T, S, V:
    R = indices[:]
    T = []
    S = np.array(values_seeding)
    V = np.array(values)

    # Extract boundary:
    B = np.ones(len(V))
    B[indices] = 2

    borders, foo1, foo2 = extract_borders(range(len(B)), B, neighbor_lists)

    #---------------------------------------------------------------------------
    # Initialize seeds with thresholded vertices:
    #---------------------------------------------------------------------------
    if use_threshold:

        # Threshold at the median depth or within maximum values in boundary:
        threshold = np.median(S[indices]) #+ np.std(S[indices])
        indices_high = [x for x in indices if S[x] >= threshold]

        # Make sure threshold is within the maximum values of the boundary:
        B = np.ones(len(S))
        B[indices] = 2
        borders, foo1, foo2 = extract_borders(range(len(B)), B, neighbor_lists)
        borders = [x for x in borders if S[x] != -1]
        if list(frozenset(indices_high).intersection(borders)):
            threshold = np.max(S[borders]) + np.std(S[borders])
            indices_high = [x for x in indices if S[x] >= threshold]

        # Extract threshold boundary vertices as seeds:
        B = -1 * np.ones(len(S))
        B[indices_high] = 2
        seeds, foo1, foo2 = extract_borders(range(len(S)), B, neighbor_lists)

    # Or initialize seeds with the maximum value point:
    else:
        seeds = [indices[np.argmax(S[indices])]]
        indices_high = []

    #---------------------------------------------------------------------------
    # Segment the mesh from the seeds iteratively toward the boundary:
    #---------------------------------------------------------------------------
    R = list(frozenset(R).difference(indices_high))
    R = list(frozenset(R).difference(seeds))
    segments = segment_rings(R, seeds, neighbor_lists, step=1)

    # Run tracks from the seeds through the segments toward the boundary:
    print('    Track through {0} vertices in {1} segments from threshold {2:0.3f}'.
          format(len(R), len(segments), threshold))
    for seed in seeds:
        track = track_to_border(seed, segments, neighbor_lists, V, borders)
        if track:
            T.append(track)

    #---------------------------------------------------------------------------
    # Filter the tracks in two ways:
    # 1. Keep tracks that have a median track value above the background value.
    # 2. Filter tracks so that there are no two track endpoint vertices
    # within a given number of edges between each other.  We select the track
    # with higher median value when their endpoints are close.
    #---------------------------------------------------------------------------
    filter_tracks = True
    if filter_tracks:

        # Compute median track values:
        Tvalues = []
        for track in T:
            Tvalues.append(np.median(V[track]))

        # Keep tracks that have a median track value above the background value:
        #background = np.median(V[[x for lst in segments for x in lst]])
        background = np.median(V[indices])
        T2 = []
        T2values = []
        for itrack, track in enumerate(T):
            if Tvalues[itrack] > background:
                T2.append(track)
                T2values.append(Tvalues[itrack])
        T = T2
        Tvalues = T2values

        # Gather endpoint vertex indices:
        E = []
        for track in T:
            E.append(track[-1])

        # Loop through endpoints:
        E2 = []
        T2 = []
        while E:

            # Find nearby endpoints to first endpoint:
            near = segment(indices, neighbor_lists, min_region_size=1,
                           seed_lists=[[E[0]]], keep_seeding=False,
                           spread_within_labels=False, labels=[],
                           label_lists=[], values=[], max_steps=min_edges)
            near = [i for i,x in enumerate(near) if x != -1]
            E_near = [x for x in E if x in near]
            if len(E_near) > 1:

                # Select endpoint with the maximum median track value:
                Inear = [i for i,x in enumerate(E) if x in E_near]
                Imax = np.argmax([Tvalues[x] for x in Inear])
                E2.append(E[Imax])
                T2.append(T[Imax])

                # Remove nearby points for the next loop:
                E = [x for i,x in enumerate(E) if i not in Inear]
                T = [x for i,x in enumerate(T) if i not in Inear]
                Tvalues = [x for i,x in enumerate(Tvalues) if i not in Inear]

            # Otherwise keep endpoint:
            else:
                E2.append(E[0])
                T2.append(T[0])
                E = E[1::]
                T = T[1::]
                Tvalues = Tvalues[1::]

        E = E2
        T = T2

    return E, T

#---------------
# Connect points
#---------------
def connect_points(indices_points, indices, L, neighbor_lists):
    """
    Connect vertices in a surface mesh to create a curve.

    The goal of this algorithm is to assign each vertex a locally optimal
    Hidden Markov Measure Field (HMMF) value and to connect vertices according
    to a cost function that penalizes vertices that do not have high likelihood
    values and have HMMF values different than their neighbors.

    We initialize the HMMF values with likelihood values normalized to the
    interval (0.5, 1.0] (to guarantee correct topology) and take those values
    that are greater than the likelihood threshold (1 for each anchor point).

    We iteratively update each HMMF value if it is near the likelihood
    threshold such that a H_step makes it cross the threshold,
    and the vertex is a "simple point" (its addition/removal alters topology).

    Parameters for computing the cost and cost gradients:

        ``wL``: weight influence of likelihood on the cost function
        ``wN``: weight influence of neighbors on the cost function
        ``H_step``: the amount that the HMMF values are incremented

    Parameters
    ----------
    indices_points : list of integers
        indices of vertices to connect (should contain > 1)
    indices : list of integers
        indices of vertices through which to connect points
    L : numpy array of floats
        likelihood values for all vertices in mesh
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex

    Returns
    -------
    skeleton : numpy array of integers
        indices to vertices

    Examples
    --------
    >>> # Connect vertices according to likelihood values in a single fold
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_vtk, read_scalars, \
    >>>                                     read_faces_points, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.features.fundi import find_endpoints, connect_points
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> # Get neighbor_lists, scalars
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> # Select a single fold:
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> faces, lines, indices, points, npoints, fold, name, input_vtk = read_vtk(fold_file)
    >>> # Test with pre-computed endpoints:
    >>> #endpoints_file = os.path.join(path, 'tests', 'connect_points_test1.vtk')
    >>> #endpoints_file = os.path.join(path, 'tests', 'connect_points_test2.vtk')
    >>> #endpoints, name = read_scalars(endpoints_file)
    >>> #indices_points = [i for i,x in enumerate(endpoints) if x > 1]
    >>> endpoints_file = os.path.join(path, 'tests', 'connect_points_test3.vtk')
    >>> endpoints, name = read_scalars(endpoints_file)
    >>> max_endpoints = max(endpoints)
    >>> indices_points = [i for i,x in enumerate(endpoints) if x == max_endpoints]
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> L, name = read_scalars(likelihood_file,True,True)
    >>> #
    >>> H = connect_points(indices_points, indices, L, neighbor_lists)
    >>> #
    >>> # View:
    >>> H[indices_points] = 1.1
    >>> rewrite_scalars(likelihood_file, 'test_connect_points.vtk', H,
    >>>                 'connected_points', fold)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_connect_points.vtk')

    """
    import numpy as np
    from mindboggle.utils.mesh import topo_test, skeletonize
    from mindboggle.labels.label import extract_borders

    # Make sure argument is a numpy array
    if not isinstance(L, np.ndarray):
        L = np.array(L)

    #-------------------------------------------------------------------------
    # Parameters:
    #-------------------------------------------------------------------------
    # Cost and cost gradient parameters:
    wN_min = 0.0  # minimum neighborhood weight
    wN_max = 2.0  # maximum neighborhood weight (trust prior more for smoother fundi)
    H_step = 0.1  # step down HMMF value

    # Parameters to speed up optimization and for termination of the algorithm:
    grad_min = 0.1  # minimum gradient factor
    grad_max = 1.0  # maximum gradient factor
    slope_exp = 2
    rate_factor = 0.9
    min_cost_change = 0.0001  # minimum change in the sum of costs
    n_tries_no_change = 3  # number of loops without sufficient change
    min_count = 50  # minimum number of iterations (to overcome initial increasing costs)
    max_count = 300  # maximum number of iterations (in case no convergence)

    # Miscellaneous parameters:
    do_skeletonize = False
    print_interval = 50

    def compute_costs(likelihoods, hmmfs, hmmfs_neighbors, numbers_of_neighbors,
                      wN, Z=[]):
        """
        Cost function for penalizing unlikely fundus curve vertices.

        This cost function penalizes vertices with low fundus likelihood values,
        and whose Hidden Markov Measure Field (HMMF) values differ from
        their neighbors:

        cost = hmmf * (1.1 - likelihood) +
               wN * sum(abs(hmmf - hmmf_neighbors)) / size_neighbors

        term 1 promotes high likelihood values
        term 2 promotes smoothness of the HMMF values

        Note: 1.1 is used instead of 1 to ensure that there is a cost
              for all points, even for those with likelihoods close to 1.

        Parameters
        ----------
        likelihoods : numpy array of floats
            likelihood values in interval [0,1]
        hmmf : numpy array of floats
            HMMF values
        hmmf_neighbors : numpy array of floats
            HMMF values of neighboring vertices for each vertex
        numbers_of_neighbors : numpy array of integers
            number of neighbors for each vertex
        wN : float
            weight influence of neighbors on cost (term 2)
        Z : numpy array of integers in {0,1} (same shape as hmmf_neighbors)
            to remove zero-padded neighborhood elements after hmmf subtraction

        Returns
        -------
        costs : numpy array of floats
            cost values

        """
        import numpy as np

        if all(numbers_of_neighbors):

            # Subtract each HMMF value from its neighbors:
            diff = abs(hmmfs - hmmfs_neighbors)

            # Remove the padding from the neighbor array:
            if np.shape(Z) == np.shape(hmmfs_neighbors):
                diff = diff * Z

            # Compute the cost for each vertex:
            costs = hmmfs * (1.1 - likelihoods) + \
                    wN * np.sum(diff,axis=0) / numbers_of_neighbors
        else:
            import sys
            sys.exit('ERROR: No HMMF neighbors to compute cost.')

        return costs

    #---------------------------------------------------------------------------
    # Initialize all Hidden Markov Measure Field (HMMF) values with
    # likelihood values (except 0) normalized to the interval (0.5, 1.0]
    # (to guarantee correct topology). Assign a 1 for each anchor point.
    # This influences surrounding vertex neighborhoods.
    # Note: 0.5 is the class boundary threshold for the HMMF values.
    H = np.zeros(len(L))
    H_new = (L + 1.000001) / 2
    H_new[L == 0.0] = 0
    H_new[H_new > 1.0] = 1
    H[H_new > 0.5] = H_new[H_new > 0.5]
    H[indices_points] = 1
    H_new = H.copy()

    # Find the HMMF values for the neighbors of each vertex:
    N = neighbor_lists
    N_sizes = np.array([len(x) for x in N])
    N_array = np.zeros((max(N_sizes[indices]), len(L)))
    for index in indices:
        N_array[0:N_sizes[index], index] = N[index]
    N_array_shape = np.shape(N_array)
    N_flat = np.ravel(N_array)
    N_flat_list = N_flat.tolist()
    H_N = np.reshape(H[N_flat_list], N_array_shape)
    ind_flat = [i for i,x in enumerate(N_flat_list) if x > 0]
    len_flat = len(N_flat_list)

    # A zero in N calls H[0], so remove zero-padded neighborhood elements:
    Z = np.zeros((max(N_sizes), len(L)))
    for index in indices:
        Z[0:N_sizes[index], index] = 1

    # Assign cost values to each vertex (for indices):
    C = np.zeros(len(L))
    C[indices] = compute_costs(L[indices], H[indices], H_N[:,indices],
                               N_sizes[indices], wN_max, Z)
    H_tests = C.copy()

    # Loop until count reaches max_count or until end_flag equals zero
    # (end_flag allows the loop to continue a few times even if no change):
    count = 0
    end_flag = 0
    wN = wN_max
    gradient_factor = grad_min
    npoints = len(indices)

    while end_flag < n_tries_no_change and count < max_count:

        # Select indices with a positive HMMF value:
        V = [indices[x] for x in np.where(H[indices] > 0.0)[0]]

        # Update neighborhood H values:
        #H_N = np.reshape(H[N_flat_list], N_array_shape)
        H_N = np.zeros(len_flat)
        H_N[ind_flat] = H[N_flat[ind_flat].tolist()]
        H_N = np.reshape(H_N, N_array_shape)

        # Compute the cost gradient for the HMMF values:
        H_decr = H - H_step
        H_decr[np.where(H_decr < 0)[0]] = 0
        C_decr = compute_costs(L[V], H_decr[V], H_N[:,V], N_sizes[V], wN, Z)
        H_tests[V] = H[V] - gradient_factor * (C[V] - C_decr)

        # For each index where the HMMF value is greater than 0:
        for index in V:
            if H[index] > 0:

                # Do not update anchor point costs:
                if index not in indices_points:

                    # Update a vertex HMMF value if it is away from the threshold:
                    update = True
                    # Or if it crosses the threshold and is a topologically
                    # "simple point" (0.5 not considered part of the fundus):
                    if H[index] > 0.5 >= H_tests[index]:
                        update, n_in = topo_test(index, H_new, N)
                    elif H[index] <= 0.5 < H_tests[index]:
                        update, n_in = topo_test(index, 1 - H_new, N)
                    if update:
                        H_new[index] = H_tests[index]
        # Update the test and cost values:
        H_tests[H_tests < 0] = 0.0
        H_tests[H_tests > 1] = 1.0

        C[V] = compute_costs(L[V], H_new[V], H_N[:,V], N_sizes[V], wN, Z)

        # Sum the cost values across all vertices and tally the number
        # of HMMF values greater than the threshold.
        # After iteration 1, compare current and previous values.
        # If the values are similar, increment end_flag:
        costs = sum(C[V].tolist())
        npoints_thr = len([x for x in H[V].tolist() if x > 0.5])

        # Terminate the loop if there are insufficient changes:
        if count > 0:
            delta_cost = (costs_previous - costs) / npoints
            delta_points = npoints_thr_previous - npoints_thr
            if delta_points == 0:
                if delta_cost < min_cost_change and count > min_count:
                    end_flag += 1
            else:
                end_flag = 0

            # Display information every n_mod iterations:
            if not np.mod(count, print_interval):
                print('      Iteration {0}: {1} crossing threshold '
                      '(wN={2:0.3f}, grad={3:0.3f}, cost={4:0.3f})'.
                      format(count, delta_points, wN, gradient_factor, delta_cost))

            # Increment the gradient factor and decrement the neighborhood factor
            # so that spacing is close in early iterations and far apart later:
            factor = (count / np.round(rate_factor*max_count))**slope_exp
            if gradient_factor < grad_max:
                gradient_factor = factor * (grad_max - grad_min) + grad_min
            if wN > wN_min:
                wN = wN_max - factor * (wN_max - wN_min)

        # Reset for next iteration:
        costs_previous = costs
        npoints_thr_previous = npoints_thr
        H = H_new

        count += 1

    print('      Updated hidden Markov measure field (HMMF) values')

    # Threshold the resulting array:
    H[H > 0.5] = 1.0
    H[H <= 0.5] = 0.0
    npoints_thr = sum(H.tolist())

    # Skeletonize:
    if do_skeletonize:
        skeleton = skeletonize(H, indices_points, N)
        print('      Removed {0} points to create one-vertex-thin skeletons'.
              format(int(npoints_thr - sum(skeleton))))
    else:
        skeleton = H

    return skeleton


# Example
if __name__ == "__main__" :

    # Extract fundus from a single fold or multiple folds:
    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    from mindboggle.features.fundi import extract_fundi

    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    mean_curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')
    likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')

    single_fold = False
    if single_fold:
        folds_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    else:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    folds, name = read_scalars(folds_file, return_first=True, return_array=True)

    fundi, n_fundi, endpoints, fundi_file = extract_fundi(folds, depth_file,
        likelihoods_or_file, save_file=True)

    # View:
    plot_vtk('fundi.vtk')


    """
    # Connect vertices according to likelihood values in a single fold
    import os
    from mindboggle.utils.io_vtk import read_vtk, read_scalars, \
                                        read_faces_points, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.features.fundi import find_endpoints, connect_points
    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    # Get neighbor_lists, scalars
    faces, points, npoints = read_faces_points(depth_file)
    neighbor_lists = find_neighbors(faces, npoints)
    # Select a single fold:
    fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    faces, lines, indices, points, npoints, fold, name, input_vtk = read_vtk(fold_file)
    # Test with pre-computed endpoints:
    #endpoints_file = os.path.join(path, 'tests', 'connect_points_test1.vtk')
    #endpoints_file = os.path.join(path, 'tests', 'connect_points_test2.vtk')
    #endpoints, name = read_scalars(endpoints_file)
    #indices_points = [i for i,x in enumerate(endpoints) if x > 1]
    endpoints_file = os.path.join(path, 'tests', 'connect_points_test3.vtk')
    endpoints, name = read_scalars(endpoints_file)
    max_endpoints = max(endpoints)
    indices_points = [i for i,x in enumerate(endpoints) if x == max_endpoints]
    likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    L, name = read_scalars(likelihood_file,True,True)
    #
    H = connect_points(indices_points, indices, L, neighbor_lists)
    #
    # View:
    H[indices_points] = 1.1
    rewrite_scalars(likelihood_file, 'test_connect_points.vtk', H,
                        'connected_points', fold)
    from mindboggle.utils.mesh import plot_vtk
    plot_vtk('test_connect_points.vtk')
    """