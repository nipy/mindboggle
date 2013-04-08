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
    >>> # Extract fundus from a single fold:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')
    >>> likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file, return_first=True, return_array=True)
    >>> #
    >>> fundi, n_fundi, endpoints, fundi_file = extract_fundi(fold, depth_file,
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
    npoints = len(folds)
    Z = np.zeros(npoints)
    fundi = -1 * np.ones(npoints)
    fundus_endpoints = []

    # For each fold region...
    unique_fold_IDs = [x for x in np.unique(folds) if x > -1]
    n_folds = len(unique_fold_IDs)
    print("Extract a fundus from each of {0} regions...".format(n_folds))
    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Region {0}:'.format(fold_ID))

            # Find fundus points
            indices_endpoints, tracks = find_endpoints(indices_fold,
                neighbor_lists, likelihoods, step=1)

            n_endpoints = len(indices_endpoints)
            if n_endpoints > 1:
                fundus_endpoints.extend(indices_endpoints)
                t2 = time()

                # Connect fundus points and extract fundi:
                print('    Connect {0} fundus points...'.format(n_endpoints))
                B = connect_points(indices_endpoints, indices_fold,
                                   likelihoods, neighbor_lists)
                indices_skeleton = [i for i,x in enumerate(B) if x > 0]
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
        rewrite_scalars(depth_file, fundi_file, [fundi],
                        ['fundi'], folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundus_endpoints, fundi_file

def segment_rings(region, seeds, neighbor_lists, step=1):
    """
    Segment a region of surface mesh iteratively toward its edges.

    Store the concentric segments for use in constructing tracks.

    Parameters
    ----------
    region : list of integers
        indices of region vertices to segment (such as a fold)
    seeds : list of integers
        indices of seed vertices
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    step : integer
        number of segmentation steps before assessing segments

    Returns
    -------
    segments : list of lists of integers
        indices to vertices for each concentric segment

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.labels.label import extract_borders
    >>> from mindboggle.features.fundi import segment_rings
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> #values, name = read_scalars(likelihood_file, True, True)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> values, name = read_scalars(depth_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> region = [i for i,x in enumerate(fold) if x != -1]
    >>> seeds = [region[np.argmax(values[region])]]
    >>> #
    >>> segments = segment_rings(region, seeds, neighbor_lists, step=1)
    >>> #
    >>> # View:
    >>> S = -1 * np.ones(len(values))
    >>> for i, segment in enumerate(segments):
    >>>     S[segment] = i
    >>> rewrite_scalars(depth_file, 'segment_rings.vtk', S, 'segment_rings', fold)
    >>> plot_vtk('segment_rings.vtk')
    >>> # Store:
    >>> #import pickle
    >>> #pickle.dump(segments, open('segments_depth_fold11.pkl', "wb" ))

    """
    from mindboggle.labels.segment import segment

    segments = []
    while seeds:

        # Segment step-wise starting from seeds and through the region:
        seeds_plus_new = segment(region, neighbor_lists, min_region_size=1,
                                 seed_lists=[seeds], keep_seeding=False,
                                 spread_within_labels=False, labels=[],
                                 label_lists=[], values=[], max_steps=step)
        seeds_plus_new = [i for i,x in enumerate(seeds_plus_new) if x != -1]

        # Store the new segment after removing the previous segment:
        region = list(frozenset(region).difference(seeds))
        seeds = list(frozenset(seeds_plus_new).difference(seeds))
        if seeds:

            # Add the new segment and remove it from the region:
            segments.append(seeds)
            region = list(frozenset(region).difference(seeds))

    return segments

def track_to_border(seed, segments, neighbor_lists, mesh_values, borders):
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
    mesh_values : numpy array of floats
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
    >>> segments_file = os.path.join(path, 'tests', 'segments_depth_fold11.pkl')
    >>> segments = pickle.load(open(segments_file, 'rb'))
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> mesh_values, name = read_scalars(depth_file, True, True)
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> region = [i for i,x in enumerate(fold) if x != -1]
    >>> seed = region[np.argmax(mesh_values[region])]
    >>> # Extract boundary:
    >>> D = np.ones(len(mesh_values))
    >>> D[region] = 2
    >>> borders, foo1, foo2 = extract_borders(range(len(mesh_values)), D, neighbor_lists)
    >>> #
    >>> track = track_to_border(seed, segments, neighbor_lists, mesh_values, borders)
    >>> #
    >>> # View:
    >>> T = -1 * np.ones(len(mesh_values))
    >>> T[track] = 1
    >>> rewrite_scalars(depth_file, 'track_to_border.vtk', T, 'track', fold)
    >>> plot_vtk('track_to_border.vtk')

    """
    import numpy as np

    track = []
    for isegment, segment in enumerate(segments):

        # Find the seed's neighborhood N in the segment:
        N = neighbor_lists[seed]
        N_segment = list(frozenset(N).intersection(segment))
        if N:

            # Add the neighborhood vertex with the maximum value to the track:
            if N_segment:
                seed = N_segment[np.argmax(mesh_values[N_segment])]
                track.append(seed)

                # If the track has run into the region's border, return the track:
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
                        if np.max(mesh_values[N_next]) > max_bridge:
                            seed = N_next[np.argmax(mesh_values[N_next])]
                            bridge = [Np, seed]
                            max_bridge = np.max(mesh_values[N_next])
                if bridge:
                    track.extend(bridge)

                    # If the track has run into the region's border, return the track:
                    if seed in borders:
                        return track

        # If there is no neighborhood for the seed, return the track:
        else:
            return track

#------------------------------------------------------------------------------
# Find endpoints
#------------------------------------------------------------------------------
def find_endpoints(indices, neighbor_lists, values, step=1):
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
    step : integer
        number of segmentation steps before assessing segments

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
    >>> values, name = read_scalars(values_file, True, True)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> step = 1
    >>> min_size = 50
    >>> #
    >>> #-----------------------------------------------------------------------
    >>> # Find endpoints on a single fold:
    >>> #-----------------------------------------------------------------------
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != -1]
    >>> indices_endpoints, tracks = find_endpoints(indices, neighbor_lists, \
    >>>                                            values, step)
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
    >>> plot_each_fold = False
    >>> folds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> fold_numbers = [x for x in np.unique(folds) if x != -1]
    >>> nfolds = len(fold_numbers)
    >>> endpoints = []
    >>> for ifold, fold_number in enumerate(fold_numbers):
    >>>     print('Fold {0} ({1} of {2})'.format(int(fold_number), ifold+1, nfolds))
    >>>     indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>>     if len(indices) > min_size:
    >>>         indices_endpoints, tracks = find_endpoints(indices, \
    >>>             neighbor_lists, values, step)
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
    from mindboggle.features.fundi import find_endpoints

    # Initialize R, T, S, V:
    R = indices[:]
    T = []
    S = np.array(values)
    V = np.copy(S)

    # Extract boundary:
    B = np.ones(len(V))
    B[indices] = 2

    borders, foo1, foo2 = extract_borders(range(len(B)), B, neighbor_lists)

    # Initialize seeds with the boundary of a thresholded region:
    use_threshold = True
    if use_threshold:
        threshold = np.max(S[borders]) + np.std(S[indices])
        indices_high = [x for x in indices if S[x] >= threshold]
        B = np.ones(len(S))
        B[indices_high] = 2
        seeds, foo1, foo2 = extract_borders(range(len(S)), B, neighbor_lists)
        R = list(frozenset(R).difference(indices_high))
    # Initialize P with the maximum value point:
    else:
        Imax = indices[np.argmax(S[indices])]
        seeds = [Imax]
        R.remove(Imax)

    # Segment the mesh iteratively toward its edges:
    segments = segment_rings(R, seeds, neighbor_lists, step=1)

    # Run recursive function track() to return track segments:
    print('  Track through {0} vertices in {1} segments from threshold {2:0.3f}'.
          format(len(R), len(segments), threshold))
    for seed in seeds:
        track = track_to_border(seed, segments, neighbor_lists, V, borders)
        if track:
            T.append(track)

    # Filter tracks by the mean mesh value:
    filter_tracks = True
    if filter_tracks:
        min_track_value = np.mean(V[indices])
        T2 = []
        E = []
        for track in T:
            track_value = np.mean(V[track])
            if track_value > min_track_value:
                T2.append(track)
                E.append(track[-1])
        T = T2

    # Extract endpoints from tracks:
    E = []
    for track in T:
        if track[-1] not in E:
            E.append(track[-1])

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

        ``H_step``: the amount that the HMMF values are H_step'd

    Parameters to speed up optimization and terminate the algorithm:

        ``min_H``: minimum HMMF value to fix very low values

        ``min_change``: minimum change in the sum of costs

        ``n_tries_no_change``: #times the loop can continue even without any change

        ``max_count``: maximum #iterations

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
    print_interval = 100

    #---------------------------------------------------------------------------
    # Cost function:
    #---------------------------------------------------------------------------
    def compute_cost(likelihood, hmmf, hmmf_neighbors, wN):
        """
        Cost function for penalizing unlikely fundus curve vertices.

        This cost function penalizes vertices with low fundus likelihood values,
        and whose Hidden Markov Measure Field (HMMF) values differ from
        their neighbors:

        cost = hmmf * (1.1 - likelihood) +
               wN * sum(abs(hmmf - hmmf_neighbors)) / len(hmmf_neighbors)

        term 1 promotes high likelihood values
        term 2 promotes smoothness of the HMMF values

        Note: 1.1 is used instead of 1 to ensure that there is a cost
              for all points, even for those with likelihoods close to 1.

        Parameters
        ----------
        likelihood : float
            likelihood value in interval [0,1]
        hmmf : float
            HMMF value
        hmmf_neighbors : numpy array of floats
            HMMF values of neighboring vertices
        wN : float
            weight influence of neighbors on cost (term 2)

        Returns
        -------
        cost : float

        """
        import numpy as np

        # Make sure arguments are numpy arrays
        if not isinstance(hmmf_neighbors, np.ndarray):
            hmmf_neighbors = np.array(hmmf_neighbors)

        if hmmf_neighbors.size:
            cost = hmmf * (1.1 - likelihood) + \
                   wN * sum(abs(hmmf - hmmf_neighbors)) / hmmf_neighbors.size
        else:
            import sys
            sys.exit('ERROR: No HMMF neighbors to compute cost.')

        return cost

    #---------------------------------------------------------------------------
    # Initialize all Hidden Markov Measure Field (HMMF) values with
    # likelihood values (except 0) normalized to the interval (0.5, 1.0]
    # (to guarantee correct topology). Assign a 1 for each anchor point.
    # This influences surrounding vertex neighborhoods.
    # Note: 0.5 is the class boundary threshold for the HMMF values.
    npoints = len(indices)
    C = np.zeros(len(L))
    H = C.copy()
    H_init = (L + 1.000001) / 2
    H_init[L == 0.0] = 0
    H_init[H_init > 1.0] = 1
    H[H_init > 0.5] = H_init[H_init > 0.5]
    #
    H[indices_points] = 1

    # Neighbors for each vertex
    N = neighbor_lists

    # Assign cost values to each vertex
    C[indices] = [compute_cost(L[i], H[i], H[N[i]], wN_max) for i in indices]

    # Loop until count reaches max_count or until end_flag equals zero
    # (end_flag is used to allow the loop to continue even if there is
    #  no change for n_tries_no_change times)
    count = 0
    end_flag = 0
    H_new = H.copy()
    wN = wN_max
    gradient_factor = grad_min
    while end_flag < n_tries_no_change and count < max_count:

        # For each index
        for index in indices:

            if H[index] > 0:

                # Do not update anchor point costs
                if index not in indices_points:

                    # Compute the cost gradient for the HMMF value
                    H_down = max([H[index] - H_step, 0])
                    cost_down = compute_cost(L[index], H_down, H[N[index]], wN)
                    H_test = H[index] - gradient_factor * (C[index] - cost_down)

                    # Update the HMMF value if near the threshold
                    # such that a step makes it cross the threshold,
                    # and the vertex is a "simple point"
                    # Note: H_new[index] is not changed yet since
                    #       topo_test() only considers its neighbors
                    if H[index] >= 0.5 >= H_test:
                        update, n_in = topo_test(index, H_new, N)
                    elif H[index] <= 0.5 <= H_test:
                        update, n_in = topo_test(index, 1 - H_new, N)

                    # Update the HMMF value if far from the threshold
                    else:
                        update = True

                    # Update the HMMF and cost values
                    if update:
                        if H_test < 0:
                            H_test = 0.0
                        elif H_test > 1:
                            H_test = 1.0
                        H_new[index] = H_test
                        C[index] = compute_cost(L[index],
                                                H_new[index], H[N[index]], wN)

        # Sum the cost values across all vertices and tally the number
        # of HMMF values greater than the threshold.
        # After iteration 1, compare current and previous values.
        # If the values are similar, increment end_flag.
        costs = sum(C)
        npoints_thr = sum([1 for x in H if x > 0.5])

        # Terminate the loop if there are insufficient changes
        if count > 0:
            delta_cost = (costs_previous - costs) / npoints
            delta_points = npoints_thr_previous - npoints_thr
            if delta_points == 0:
                if delta_cost < min_cost_change and count > min_count:
                    end_flag += 1
            else:
                end_flag = 0

            # Display information every n_mod iterations
            if not np.mod(count, print_interval):
                print('      Iteration {0}: {1} points crossing threshold (wN={2:0.3f}, grad={3:0.3f}, cost={4:0.3f})'.
                      format(count, delta_points, wN, gradient_factor, delta_cost))

            # Increment the gradient factor and
            # decrement the neighborhood factor
            # so that the spacing is close in early iterations
            # and far apart in later increments.
            factor = (count / np.round(rate_factor*max_count))**slope_exp
            if gradient_factor < grad_max:
                gradient_factor = factor * (grad_max - grad_min) + grad_min
            if wN > wN_min:
                wN = wN_max - factor * (wN_max - wN_min)

        # Reset for next iteration
        costs_previous = costs
        npoints_thr_previous = npoints_thr
        H = H_new

        count += 1

    print('      Updated hidden Markov measure field (HMMF) values')

    # Threshold the resulting array
    H[H > 0.5] = 1
    H[H <= 0.5] = 0
    npoints_thr = sum(H)

    # Skeletonize
    if do_skeletonize:
        skeleton = skeletonize(H, indices_points, N)
        print('      Removed {0} points to create one-vertex-thin skeletons'.
              format(int(npoints_thr - sum(skeleton))))
    else:
        skeleton = H

    return skeleton


# Example
if __name__ == "__main__" :

    # Extract fundus from a single fold:
    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    from mindboggle.features.fundi import extract_fundi
    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    mean_curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')
    likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    fold, name = read_scalars(fold_file, return_first=True, return_array=True)

    fundi, n_fundi, endpoints, fundi_file = extract_fundi(fold, depth_file,
        likelihoods_or_file, save_file=True)

    # View:
    plot_vtk('fundi.vtk')
