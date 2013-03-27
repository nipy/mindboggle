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
def extract_fundi(folds_or_file, depth_rescaled_file,
                  min_curvature_vector_file, likelihoods_or_file,
                  min_distance=5, thr=0.5, save_file=False):
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
    depth_rescaled_file :  string
        surface mesh file in VTK format with scalar rescaled depth values
    min_curvature_vector_file : string
        surface mesh file in VTK format with scalar values
    likelihoods_or_file : list or string
        fundus likelihood values or name of VTK file with the scalar values
    min_fold_size : int
        minimum fold size (number of vertices)
    min_distance :  int
        minimum distance between "anchor points"
    thr :  float
        likelihood threshold
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundi : list of integers
        fundus numbers for all vertices (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundus_anchor_indices : list of integers
        indices to fundus anchor vertices
    fundi_file : string (if save_file)
        name of output VTK file with fundus numbers (-1 for non-fundus vertices)
        and likelihood values for sulcus vertices (separate scalars)

    Examples
    --------
    >>> # Extract fundus from a single fold:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_rescaled_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')
    >>> likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> # Select a single fold
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file, return_first=True, return_array=True)
    >>>
    >>> fundi, n_fundi, likelihoods, fundus_anchor_indices, fundi_file = extract_fundi(fold,
    >>>     depth_rescaled_file, min_curv_vec_file,
    >>>     likelihoods_or_file, min_distance=5, thr=0.5, save_file=True)
    >>>
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('fundi.vtk')

    """
    import os
    import sys
    import numpy as np
    from time import time

    from mindboggle.features.fundi import find_endpoints, connect_points
    from mindboggle.features.likelihood import compute_likelihood
    from mindboggle.utils.mesh import find_neighbors, skeletonize, extract_skeleton_endpoints
    from mindboggle.utils.io_vtk import read_scalars, read_faces_points, \
        read_vtk, rewrite_scalars

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

    faces, points, npoints = read_faces_points(depth_rescaled_file)
    min_directions = np.loadtxt(min_curvature_vector_file)

    # Initialize variables:
    t1 = time()
    count = 0
    neighbor_lists = find_neighbors(faces, npoints)
    points = np.array(points)
    Z = np.zeros(npoints)
    fundi = -1 * np.ones(npoints)
    fundus_anchor_indices = []
    if compute_likelihoods:
        likelihoods = np.copy(fundi)

    # For each fold region...
    unique_fold_IDs = [x for x in np.unique(folds) if x > -1]
    n_folds = len(unique_fold_IDs)
    print("Extract a fundus from each of {0} regions...".format(n_folds))
    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:

            print('  Region {0}:'.format(fold_ID))

            fold_likelihoods = likelihoods[indices_fold]

            # Find fundus points
            ########indices, neighbor_lists, likelihoods, step_size=5
            fold_indices_endpoints = find_endpoints(points[indices_fold],
                                                    fold_likelihoods,
                                                    min_directions[indices_fold],
                                                    min_distance, thr)
            indices_endpoints = [indices_fold[x] for x in fold_indices_endpoints]
            n_endpoints = len(indices_endpoints)
            if n_endpoints > 1:
                t2 = time()
                likelihoods_fold = Z.copy()
                likelihoods_fold[indices_fold] = fold_likelihoods

                # Connect fundus points and extract fundi:
                print('    Connect {0} fundus points...'.format(n_endpoints))
                B = connect_points(indices_endpoints, indices_fold,
                                   likelihoods_fold, neighbor_lists)
                indices_skeleton = [i for i,x in enumerate(B) if x > 0]
                fundus_anchor_indices.extend(indices_endpoints)

                if len(indices_skeleton) > 1:
                    fundi[indices_skeleton] = fold_ID
                    count += 1
                    print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                          format(n_endpoints, time() - t2))

    n_fundi = count
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.format(n_fundi, time() - t1))

    #---------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name
    #---------------------------------------------------------------------------
    fundi = fundi.tolist()

    if save_file:
        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(depth_rescaled_file, fundi_file, [fundi],
                        ['fundi'], folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundus_anchor_indices, fundi_file


#--------------
# Cost function
#--------------
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

    Note: 1.1 is used instead of 1 to ensure that there is a cost for all points,
          even for those with likelihoods close to 1.

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

#------------------------------------------------------------------------------
# Find special points that are not too close together
#------------------------------------------------------------------------------
def find_endpoints(indices, neighbor_lists, likelihoods, step_size=5):
    """
    Find endpoints in a region of connected vertices.

    These points are intended to serve as endpoints of fundus curves
    running along high-likelihood paths.  This algorithm iteratively
    propagates from an initial high-likelihood seed toward the boundary
    of the region, selecting the highest likelihood boundary point
    from each terminal segment as an endpoint.

    Steps ::

        Initialize:
            R: Region/remaining vertices to segment (initially fold vertices)
            P: Previous segment (initially the maximum likelihood point)
            N: New/next segment (initially P)
            E: indices to endpoint vertices (initially empty)
        For each segment P, run recursive function creep():
            Propagate P into R, and call these new vertices N.
            If N is empty:
                Choose highest likelihood point in P as endpoint.
                Return endpoints E and remaining vertices R.
            else:
                Identify N_i different segments of N.
                For each segment N_i:
                    If N_i large enough or if max(i)==1:
                        Call creep() with new arguments.
                Return endpoints E and remaining vertices R.

    Parameters
    ----------
    indices : list of integers
        indices of the vertices to segment (such as a fold in a surface mesh)
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    likelihoods : numpy array of floats
        fundus likelihood values for all vertices
    step_size : integer
        number of segmentation steps before assessing segments

    Returns
    -------
    indices_endpoints : list of integers
        indices of surface mesh vertices that are endpoints

    Examples
    --------
    >>> # Find endpoints on a single fold:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.fundi import find_endpoints
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> likelihoods, name = read_scalars(likelihood_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(likelihood_file)
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != -1]
    >>> step_size = 1
    >>> #
    >>> indices_endpoints = find_endpoints(indices, neighbor_lists, likelihoods, step_size)
    >>> #
    >>> # Write results to VTK file and view:
    >>> likelihoods[indices_endpoints] = max(likelihoods) + 0.1
    >>> rewrite_scalars(likelihood_file, 'find_endpoints.vtk',
    >>>                 likelihoods, 'endpoints_on_likelihoods_in_folds', fold)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('find_endpoints.vtk')

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if isinstance(likelihoods, list):
        likelihoods = np.array(likelihoods)

    # Recursive function for segmenting and finding endpoints:
    def creep(R, P, N, E, L, step_size, neighbor_lists):
        """
        Recursively segment a mesh, creeping toward its edges to find endpoints.

        Steps ::

            Propagate P into R, and call these new vertices N.
            If N is empty:
                Choose highest likelihood point in P as endpoint.
                Return endpoints E and remaining vertices R.
            else:
                Identify N_i different segments of N.
                For each segment N_i:
                    If N_i large enough or if max(i)==1:
                        Call creep() with new arguments.
                Return endpoints E and remaining vertices R.

        Parameters
        ----------
        R : list of integers
            indices of vertices to segment (such as a fold)
        P : list of integers
            indices of previous segment vertices
        N : list of integers
            indices of new/next segment vertices
        E: list of integers
            indices to endpoint vertices
        L : numpy array of floats
            likelihood values for all vertices
        step_size : integer
            number of segmentation steps before assessing segments
        neighbor_lists : list of lists of integers
            indices to neighboring vertices for each vertex

        Returns
        -------
        R : list of integers
            remaining vertices to segment
        P : list of integers
            previous segment vertices
        N : list of integers
            new segment vertices
        E: list of integers
            endpoint vertices

        """
        import numpy as np

        from mindboggle.labels.segment import segment

        min_size = 3

        #-----------------------------------------------------------------------
        # Propagate P into R, and call these new vertices N:
        #-----------------------------------------------------------------------
        R_segments = segment(R, neighbor_lists, min_region_size=1,
                             seed_lists=[P], keep_seeding=False, spread_within_labels=False,
                             labels=[], label_lists=[], values=[], max_steps=step_size)
        N = [i for i,x in enumerate(R_segments) if x != -1]

        # Remove P (seeds) from N, and remove P and N from R:
        N = list(frozenset(N).difference(P))
        R = list(frozenset(R).difference(P))
        R = list(frozenset(R).difference(N))

        #-----------------------------------------------------------------------
        # If N is empty, return endpoints:
        #-----------------------------------------------------------------------
        if not N:

            # Choose highest likelihood point in P as endpoint:
            E.append(P[np.argmax(L[P])])

            # Return endpoints E and remaining vertices R:
            return R, P, N, E

        #-----------------------------------------------------------------------
        # If N is not empty, continue segmenting recursively:
        #-----------------------------------------------------------------------
        else:

            # Identify N_i different segments of N:
            N_segments = segment(N, neighbor_lists, min_region_size=1)
            unique_N = [x for x in np.unique(N_segments) if x!=-1]
            n_segments = len(unique_N)

            # For each segment N_i:
            for n in unique_N:
                N_i = [i for i,x in enumerate(N_segments) if x==n]

                # If N_i large enough or if max(i)==1:
                if len(N_i) >= min_size or n_segments==1:

                    # Call creep() with new arguments:
                    R, P, N, E = creep(R, N_i, N_i, E, L, step_size, neighbor_lists)

            # Return endpoints E and remaining vertices R:
            return R, P, N, E


    # Initialize old and new segments P and N with the maximum likelihood point:
    maxL = indices[np.argmax(likelihoods[indices])]

    # For each threshold:
    indices_endpoints = []
    thresholds = [0.5]
    for thr in thresholds:

        # Run recursive function creep() to return endpoints:
        R = indices
        P = [maxL]
        N = [maxL]
        E = indices_endpoints
        L = likelihoods
        R, P, N, E = creep(R, P, N, E, L, step_size, neighbor_lists)

        indices_endpoints.extend(E)

    return indices_endpoints

#---------------
# Connect points
#---------------
def connect_points(indices_endpoints, indices, L, neighbor_lists):
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
    indices_endpoints : list of integers
        indices of vertices to connect (should contain > 1)
    indices : list of integers
        indices of vertices
    L : numpy array of floats
        likelihood values: #vertices in mesh x 1
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
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.features.fundi import find_endpoints, connect_points
    >>> from mindboggle.features.likelihood import compute_likelihood
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_rescaled_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> # Get neighbor_lists, scalars
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_rescaled_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> # Select a single fold:
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file, return_first=True, return_array=True)
    >>> # Test with pre-computed endpoints:
    >>> endpoints_file = os.path.join(path, 'tests', 'connect_points_test1.vtk')
    >>> #endpoints_file = os.path.join(path, 'tests', 'connect_points_test2.vtk')
    >>> endpoints, name = read_scalars(endpoints_file)
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> L, name = read_scalars(likelihood_file,True,True)
    >>> indices_endpoints =[i for i,x in enumerate(endpoints) if x>1 if folds[i]==fold_ID]
    >>> #
    >>> H = connect_points(indices_endpoints, indices, L, neighbor_lists)
    >>> #
    >>> # View:
    >>> H[indices_endpoints] = 1.1
    >>> rewrite_scalars(depth_rescaled_file, 'test_connect_points.vtk', H,
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
    # Parameters
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
    #-------------------------------------------------------------------------

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
    H[indices_endpoints] = 1
    print_interval = 100

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
                if index not in indices_endpoints:

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
                print('      Iteration {0}: {1} points crossing threshold (wN={2}, grad={3}, cost={4})'.
                      format(count, delta_points, wN, gradient_factor, delta_cost))
                #print('      H[0]: {0}'.format(H[indices]))

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
        skeleton = skeletonize(H, indices_endpoints, N)
        print('      Removed {0} points to create one-vertex-thin skeletons'.
              format(int(npoints_thr - sum(skeleton))))
    else:
        skeleton = H

    return skeleton



# Example
if __name__ == "__main__" :

    # Extract fundus from a single fold:
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.features.fundi import extract_fundi
    path = os.environ['MINDBOGGLE_DATA']
    depth_rescaled_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    min_curv_vec_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.dir.txt')

    # Select a single fold
    folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    folds, name = read_scalars(folds_file)
    fold_ID = 3
    indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    fold_array = -1 * np.ones(len(folds))
    fold_array[indices_fold] = 1

    fundi, n_fundi, likelihoods, fundus_anchor_indices, fundi_file = extract_fundi(fold_array,
        depth_rescaled_file, min_curv_vec_file,
        likelihoods_or_file, min_distance=5, thr=0.5, save_file=True)

    # View:
    from mindboggle.utils.mesh import plot_vtk
    #plot_vtk('fundi.vtk')

    # View overlay of likelihood values, fundus, and anchor points:
    overlay = np.array(likelihoods_or_file)
    Ifundi = [i for i,x in enumerate(fundi) if x > -1]
    overlay[Ifundi] = 0.25
    overlay[fundus_anchor_indices] = 0
    rewrite_scalars(depth_rescaled_file, 'overlay.vtk', overlay,
                    'likelihoods_fundi_endpoints', fold_array)
    plot_vtk('overlay.vtk')
