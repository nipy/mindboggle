#!/usr/bin/env python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
<<<<<<< HEAD
Arno Klein, 2013  .  arno@mindboggle.info  .  www.binarybottle.com
=======
Yrjo Hame, 2012-2013  .  yrjo.hame@gmail.com
Arno Klein, 2012-2013  .  arno@mindboggle.info  .  www.binarybottle.com
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

<<<<<<< HEAD

def extract_fundi(folds, curv_file, depth_file, min_separation=10,
                  erode_ratio=0.1, erode_min_size=1, save_file=False):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a fold.

    Steps ::
        1. Find fundus endpoints (outer anchors) with find_outer_anchors().
        2. Include inner anchor points.
        3. Connect anchor points using connect_points_erosion();
           inner anchors are removed if they result in endpoints.

    Parameters
    ----------
    folds : numpy array or list of integers
        fold number for each vertex
    curv_file :  string
        surface mesh file in VTK format with mean curvature values
    depth_file :  string
        surface mesh file in VTK format with rescaled depth values
    likelihoods : list of integers
        fundus likelihood value for each vertex
    min_separation : integer
        minimum number of edges between inner/outer anchor points
    erode_ratio : float
        fraction of indices to test for removal at each iteration
        in connect_points_erosion()
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundus_per_fold : list of integers
        fundus numbers for all vertices, labeled by fold
        (-1 for non-fundus vertices)
    n_fundi_in_folds :  integer
        number of fundi
    fundus_per_fold_file : string (if save_file)
        output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> single_fold = True
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.features.fundi import extract_fundi
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'travel_depth_rescaled.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> if single_fold:
    >>>     fold_number = 2 #11
    >>>     folds[folds != fold_number] = -1
    >>> min_separation = 10
    >>> erode_ratio = 0.10
    >>> erode_min_size = 10
    >>> save_file = True
    >>> o1, o2, fundus_per_fold_file = extract_fundi(folds, curv_file,
    ...     depth_file, min_separation, erode_ratio, erode_min_size, save_file)
    >>> #
    >>> # View:
    >>> plot_surfaces(fundi_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np
    from time import time

    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.compute import median_abs_dev
    from mindboggle.utils.paths import find_max_values
    from mindboggle.utils.mesh import find_neighbors_from_file, find_complete_faces
    from mindboggle.utils.paths import find_outer_anchors, connect_points_erosion

    if isinstance(folds, list):
        folds = np.array(folds)

    # Load values, inner anchor threshold, and neighbors:
    faces, u1,u2, points, npoints, curvs, u3,u4 = read_vtk(curv_file, True,True)
    depths, name = read_scalars(depth_file, True, True)
    values = curvs * depths
    values0 = [x for x in values if x > 0]
    thr = np.median(values0) + 2 * median_abs_dev(values0)
    neighbor_lists = find_neighbors_from_file(curv_file)

    #-------------------------------------------------------------------------
    # Loop through folds:
    #-------------------------------------------------------------------------
    t1 = time()
    skeletons = []
    unique_fold_IDs = [x for x in np.unique(folds) if x != -1]

    if len(unique_fold_IDs) == 1:
        print("Extract a fundus from 1 fold...")
    else:
        print("Extract a fundus from each of {0} folds...".
              format(len(unique_fold_IDs)))

    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Fold {0}:'.format(int(fold_ID)))

            #-----------------------------------------------------------------
            # Find outer anchor points on the boundary of the surface region,
            # to serve as fundus endpoints:
            #-----------------------------------------------------------------
            outer_anchors, tracks = find_outer_anchors(indices_fold,
                neighbor_lists, values, depths, min_separation)

            #-----------------------------------------------------------------
            # Find inner anchor points:
            #-----------------------------------------------------------------
            inner_anchors = find_max_values(points, values, min_separation, thr)

            #-----------------------------------------------------------------
            # Connect anchor points to create skeleton:
            #-----------------------------------------------------------------
            B = -1 * np.ones(npoints)
            B[indices_fold] = 1
            skeleton = connect_points_erosion(B, neighbor_lists,
                outer_anchors, inner_anchors, values,
                erode_ratio, erode_min_size, save_steps=[], save_vtk='')
            if skeleton:
                skeletons.extend(skeleton)

            #-----------------------------------------------------------------
            # Remove fundus vertices if they complete triangle faces:
            #-----------------------------------------------------------------
            Iremove = find_complete_faces(skeletons, faces)
            if Iremove:
                skeletons = list(frozenset(skeletons).difference(Iremove))

    indices = [x for x in skeletons if folds[x] != -1]
    fundus_per_fold = -1 * np.ones(npoints)
    fundus_per_fold[indices] = folds[indices]
    n_fundi_in_folds = len([x for x in np.unique(fundus_per_fold)
                             if x != -1])
    if n_fundi_in_folds == 1:
        sdum = 'fold fundus'
    else:
        sdum = 'fold fundi'
    print('  ...Extracted {0} {1}; {2} total ({3:.2f} seconds)'.
          format(n_fundi_in_folds, sdum, n_fundi_in_folds, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    if n_fundi_in_folds > 0:
        fundus_per_fold = [int(x) for x in fundus_per_fold]
        if save_file:
            fundus_per_fold_file = os.path.join(os.getcwd(),
                                                'fundus_per_fold.vtk')
            rewrite_scalars(curv_file, fundus_per_fold_file, fundus_per_fold,
                            'fundi', folds)
            if not os.path.exists(fundus_per_fold_file):
                raise(IOError(fundus_per_fold_file + " not found"))
        else:
            fundus_per_fold_file = None

    return fundus_per_fold,  n_fundi_in_folds, fundus_per_fold_file


def segment_fundi(fundus_per_fold, sulci=[], vtk_file='', save_file=False):
    """
    Segment fundi by sulcus definitions.

    Parameters
    ----------
    fundus_per_fold : list of integers
        fundus numbers for all vertices, labeled by fold
        (-1 for non-fundus vertices)
    sulci : numpy array or list of integers
        sulcus number for each vertex, used to filter and label fundi
    vtk_file : string (if save_file)
        VTK file with sulcus number for each vertex
=======
#-----------------
# Sigmoid function
#-----------------
def sigmoid(values, gain, shift):
    """
    Map values with sigmoid function to range [0,1].

    Y(t) = 1/(1 + exp(-gain*(values - shift))
    """
    import numpy as np

    # Make sure argument is a numpy array
    if type(values) != np.ndarray:
        values = np.array(values)

    return 1.0 / (1.0 + np.exp(-gain * (values - shift)))

#--------------------
# Likelihood function
#--------------------
def compute_likelihood(depths, curvatures):
    """
    Compute (fundus) curve likelihood values on a surface mesh.

    The *slope_factor* is used to compute the "gain" of the slope of sigmoidal values.

    Parameters
    ----------
    depths : numpy array of floats
        depth values in [0,1] for all vertices
    curvatures : numpy array of floats
        mean curvature values in [-1,1] for all vertives

    Returns
    -------
    likelihoods : numpy array of floats
        likelihood values for all vertices

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.features.fundi import compute_likelihood
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.avg.vtk')
    >>> depths, name = read_scalars(depth_file, return_first=True, return_array=True)
    >>> mean_curvs, name = read_scalars(mean_curv_file, return_first=True, return_array=True)
    >>>
    >>> L = compute_likelihood(depths, mean_curvs)
    >>>
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_compute_likelihood.vtk', L, 'likelihoods')
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_compute_likelihood.vtk')

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(depths, np.ndarray):
        depths = np.array(depths)
    if not isinstance(curvatures, np.ndarray):
        curvatures = np.array(curvatures)

    #==========================================
    # Normalize depth values to interval [0,1]
    # Curvature values retain their values
    # Compute the means and std. deviations
    #=========================================
    depths_norm = depths / max(depths)
    depth_avg = np.mean(depths_norm)
    depth_std = np.std(depths_norm)
    curve_avg = np.mean(curvatures)
    curve_std = np.std(curvatures)

    #==========================================
    # Find slope for depth and curvature values
    #==========================================
    # Factor influencing "gain" or "sharpness" of the sigmoidal function below
    # slope_factor = abs(np.log((1. / x) - 1))  # 2.197224577 for x = 0.9
    # gain_depth = slope_factor / (2 * depth_std)
    gain_depth = 1 / depth_std
    gain_curve = 1 / curve_std
    shift_depth = depth_avg - depth_std
    shift_curve = curve_avg

    #==========================
    # Compute likelihood values
    #==========================
    # Map values with sigmoid function to range [0,1]
    depth_sigmoid = sigmoid(depths_norm, gain_depth, shift_depth)
    curve_sigmoid = sigmoid(curvatures, gain_curve, shift_curve)

    likelihoods = depth_sigmoid * curve_sigmoid

    # Plot the sigmoid curves (does not include value distributions)
    plot_result = False
    if plot_result:
        from matplotlib import pyplot
        xdepth = np.sort(depths_norm)
        xcurve = np.sort(curvatures)
        depth_sigmoid_sort = sigmoid(xdepth, gain_depth, shift_depth)
        curve_sigmoid_sort = sigmoid(xcurve, gain_curve, shift_curve)
        sigmoids = depth_sigmoid_sort * curve_sigmoid_sort
        pyplot.plot(xdepth, depth_sigmoid_sort, 'k')
        pyplot.plot(xcurve, curve_sigmoid_sort, 'b')
        pyplot.plot(xdepth, sigmoids, 'r')
        pyplot.title('Depths, curves: (gains={0:.2f},{1:.2f}; shifts={2:.2f},{3:.2f})'.
               format(gain_depth, gain_curve, shift_depth, shift_curve))
        pyplot.show()

    return likelihoods

#--------------
# Cost function
#--------------
def compute_cost(likelihood, hmmf, hmmf_neighbors, wN):
    """
    Cost function for penalizing unlikely curve (fundus) vertices.

    This cost function penalizes vertices with low likelihood values,
    and whose Hidden Markov Measure Field (HMMF) values differ from
    their neighbors:

    cost = hmmf * (1 - likelihood) +
           wN * sum(abs(hmmf - hmmf_neighbors)) / len(hmmf_neighbors)

    term 1 promotes high likelihood values
    term 2 promotes smoothness of the HMMF values

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
    cost : ``float``

    """
    import numpy as np

    if hmmf_neighbors.size:

        # Make sure arguments are numpy arrays
        if not isinstance(hmmf_neighbors, np.ndarray):
            hmmf_neighbors = np.array(hmmf_neighbors)

        cost = hmmf * (1 - likelihood) + \
               wN * sum(abs(hmmf - hmmf_neighbors)) / len(hmmf_neighbors)
    else:
        exit('ERROR: No HMMF neighbors to compute cost.')

    return cost

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
    L : list (or array) of integers
        fundus likelihood values for all vertices (default -1)
    min_directions : numpy array of floats
        minimum directions for all vertices
    min_distance : integer
        minimum distance
    thr : float
        likelihood threshold in [0,1]

    Returns
    -------
    anchors : list of subset of surface mesh vertex indices

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.features.fundi import find_anchors
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> min_curvature_vector_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.min.dir.txt')
    >>> faces, lines, indices, points, npoints, values, name = read_vtk(depth_file)
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> min_distance = 5
    >>> thr = 0.5
    >>> #
    >>> anchors = find_anchors(points, values, min_directions, min_distance, thr)
    >>> #
    >>> # Write results to vtk file and view:
    >>> IDs = -1 * np.ones(len(min_directions))
    >>> IDs[anchors] = 1
    >>> rewrite_scalars(depth_file, 'test_find_anchors.vtk', IDs, 'anchors', IDs)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_find_anchors.vtk')

    """
    import numpy as np
    from operator import itemgetter

    # Make sure arguments are numpy arrays
    if not isinstance(points, np.ndarray):
        points = np.array(points)
    if not isinstance(min_directions, np.ndarray):
        min_directions = np.array(min_directions)

    max_distance = 2 * min_distance

    # Sort likelihood values and find indices for values above the threshold
    L_table = [[i,x] for i,x in enumerate(L)]
    L_table_sort = np.transpose(sorted(L_table, key=itemgetter(1)))[:, ::-1]
    IL = [int(L_table_sort[0,i]) for i,x in enumerate(L_table_sort[1,:])
          if x > thr]

    # Initialize anchors list with the index of the maximum likelihood value,
    # remove this value, and loop through the remaining high likelihoods
    if IL:
        anchors = [IL.pop(0)]
        for imax in IL:

            # Determine if there are any anchor points
            # near to the current maximum likelihood vertex
            i = 0
            found = 0
            while i < len(anchors) and found == 0:

                # Compute Euclidean distance between points
                D = np.linalg.norm(points[anchors[i]] - points[imax])

                # If distance less than threshold, consider the point found
                if D < min_distance:
                    found = 1
                # Compute directional distance between points if they are close
                elif D < max_distance:
                    dirV = np.dot(points[anchors[i]] - points[imax],
                                  min_directions[anchors[i]])
                    # If distance less than threshold, consider the point found
                    if np.linalg.norm(dirV) < min_distance:
                        found = 1

                i += 1

            # If there are no nearby anchor points,
            # assign the maximum likelihood vertex as an anchor point
            if not found:
                anchors.append(imax)
    else:
        anchors = []

    return anchors

#---------------
# Connect points
#---------------
def connect_points(anchors, indices, L, neighbor_lists):
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
    anchors : list of integers
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
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.features.fundi import find_anchors, connect_points, compute_likelihood
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.avg.vtk')
    >>> min_curvature_vector_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.min.dir.txt')
    >>> # Get neighbor_lists, scalars
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> points = np.array(points)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> mean_curvatures, name = read_scalars(mean_curvature_file, True, True)
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> # Select a single fold
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> fold_ID = 10
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>> fold_array = -1 * np.ones(npoints)
    >>> fold_array[indices_fold] = 1
    >>> fold_likelihoods = compute_likelihood(depths[indices_fold],
    >>>                                       mean_curvatures[indices_fold])
    >>> fold_indices_anchors = find_anchors(points[indices_fold],
    >>>     fold_likelihoods, min_directions[indices_fold], 5, 0.5)
    >>> indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
    >>> likelihoods_fold = np.zeros(len(points))
    >>> likelihoods_fold[indices_fold] = fold_likelihoods
    >>>
    >>> H = connect_points(indices_anchors, indices_fold, likelihoods_fold, neighbor_lists)
    >>>
    >>> # Write results to vtk file and view:
    >>> H[indices_anchors] = 2
    >>> rewrite_scalars(depth_file, 'test_connect_points.vtk', H, 'connected_points', fold_array)
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
    # Cost and cost gradient parameters
    wN_max = 0.5  # maximum neighborhood weight
    wN_min = 0.1  # minimum neighborhood weight
    H_step = 0.1  # step down HMMF value
    H_min = 0.5 - H_step  # minimum HMMF value to be processed

    # Parameters to speed up optimization and for termination of the algorithm
    grad_min = 0.1  # minimum gradient factor
    grad_max = 0.2  # maximum gradient factor
    min_cost_change = 0.0001  # minimum change in the sum of costs
    n_tries_no_change = 3  # number of loops without sufficient change
    max_count = 1000  # maximum number of iterations (in case no convergence)
    #-------------------------------------------------------------------------

    # Initialize all Hidden Markov Measure Field (HMMF) values with
    # likelihood values (except 0) normalized to the interval (0.5, 1.0]
    # (to guarantee correct topology). Assign a 1 for each anchor point.
    # Note: 0.5 is the class boundary threshold for the HMMF values.
    npoints = len(indices)
    C = np.zeros(len(L))
    H = C.copy()
    H_init = (L + 1.000001) / 2
    H_init[L == 0.0] = 0
    H_init[H_init > 1.0] = 1
    H[H_init > 0.5] = H_init[H_init > 0.5]
    H[anchors] = 1

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
    gradient_factor = grad_max
    while end_flag < n_tries_no_change and count < max_count:

        # For each index
        for index in indices:

            # Do not update anchor point costs
            if index not in anchors:

                # Continue if the HMMF value is greater than a minimum value
                # (fix at very low values to speed up optimization)
                if H[index] > H_min:

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
        n_points = sum([1 for x in H if x > 0.5])

        # Terminate the loop if there are insufficient changes
        if count > 0:
            delta_cost = (costs_previous - costs) / npoints
            delta_points = n_points_previous - n_points
            if delta_points == 0:
                if delta_cost < min_cost_change:
                    end_flag += 1
            else:
                end_flag = 0

            # Display information every n_mod iterations
            if not np.mod(count, 10):
                print('      Iteration {0}: {1} points crossing threshold'.
                      format(count, delta_points)) #delta_cost

            # Increment next gradient factor and decrement next neighborhood
            # weight as a function of net number of points crossing threshold
            n_change = max([abs(delta_points), 20])
            gradient_factor = grad_min + (grad_max - grad_min) / n_change
            wN = wN_max - (wN_max - wN_min) / n_change

        # Reset for next iteration
        costs_previous = costs
        n_points_previous = n_points
        H = H_new

        count += 1

    print('      Updated hidden Markov measure field (HMMF) values')

    # Threshold the resulting array
    H[H > 0.5] = 1
    H[H <= 0.5] = 0
    n_points = sum(H)

    # Skeletonize
    skeleton = skeletonize(H, anchors, N)
    print('      Removed {0} points to create one-vertex-thin skeletons'.
          format(int(n_points - sum(skeleton))))

    return skeleton

#==================
# Extract all fundi
#==================
def extract_fundi(folds, neighbor_lists, depth_file,
                  mean_curvature_file, min_curvature_vector_file,
                  min_distance=5, thr=0.5, use_only_endpoints=True,
                  compute_local_depth=True, save_file=False):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Parameters
    ----------
    folds : list or array of integers
        fold IDs (default = -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    depth_file : string
        surface mesh file in VTK format with faces and scalar values
    mean_curvature_file : string
        surface mesh file in VTK format with scalar values
    min_curvature_vector_file : string
        surface mesh file in VTK format with scalar values
    min_fold_size : int
        minimum fold size (number of vertices)
    min_distance :  int
        minimum distance between "anchor points"
    thr :  float
        likelihood threshold
    use_only_endpoints : Boolean
        use only endpoints to construct fundi (or all anchor points)?
    compute_local_depth : Boolean
        normalize depth for each fold?
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    save_file : Boolean
        save output VTK file?

    Returns
    -------
<<<<<<< HEAD
    fundus_per_sulcus : list of integers
        fundus numbers for all vertices, labeled by sulcus
        (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundus_per_sulcus_file : string (if save_file)
        output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more sulci:
    >>> single_fold = True
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.features.fundi import extract_fundi, segment_fundi
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(vtk_file, True, True)
    >>> curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'travel_depth_rescaled.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> if single_fold:
    >>>     fold_number = 2 #11
    >>>     folds[folds != fold_number] = -1
    >>> min_separation = 10
    >>> erode_ratio = 0.10
    >>> erode_min_size = 10
    >>> save_file = True
    >>> fundus_per_fold, o1, o2 = extract_fundi(folds, curv_file, depth_file, min_separation, erode_ratio, erode_min_size, save_file)
    >>> o1, o2, fundus_per_sulcus_file = segment_fundi(fundus_per_fold, sulci, vtk_file, save_file)
    >>> #
    >>> # View:
    >>> plot_surfaces(fundus_per_sulcus_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import rewrite_scalars

    if isinstance(sulci, list):
        sulci = np.array(sulci)

    #-------------------------------------------------------------------------
    # Create fundi by segmenting fold fundi with overlapping sulcus labels:
    #-------------------------------------------------------------------------
    indices = [i for i,x in enumerate(fundus_per_fold) if x != -1]
    if indices and np.size(sulci):
        fundus_per_sulcus = -1 * np.ones(len(sulci))
        fundus_per_sulcus[indices] = sulci[indices]
        n_fundi = len([x for x in np.unique(fundus_per_sulcus) if x != -1])
    else:
        fundus_per_sulcus = []
        n_fundi = 0

    if n_fundi == 1:
        sdum = 'sulcus fundus'
    else:
        sdum = 'sulcus fundi'
    print('  Segmented {0} {1}'.format(n_fundi, sdum))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    fundus_per_sulcus_file = None
    if n_fundi > 0:
        fundus_per_sulcus = [int(x) for x in fundus_per_sulcus]
        if save_file and os.path.exists(vtk_file):
            fundus_per_sulcus_file = os.path.join(os.getcwd(),
                                                  'fundus_per_sulcus.vtk')
            rewrite_scalars(vtk_file, fundus_per_sulcus_file,
                            fundus_per_sulcus, 'fundus_per_sulcus',
                            fundus_per_sulcus)
            if not os.path.exists(fundus_per_sulcus_file):
                raise(IOError(fundus_per_sulcus_file + " not found"))

    return fundus_per_sulcus, n_fundi, fundus_per_sulcus_file
=======
    fundi : array of integers
        fundus numbers for all vertices (-1 for non-fundus vertices)
    n_fundi :  int
        number of fundi
    likelihoods : array of floats
        fundus likelihood values for all vertices (zero outside folds)
    fundi_file : string (if save_file)
        name of output VTK file with fundus numbers (-1 for non-fundus vertices)
        and likelihood values for sulcus vertices (separate scalars)

    Examples
    --------
    >>> # Extract fundus from a single fold:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_faces_points, read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.avg.vtk')
    >>> min_curv_vec_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.min.dir.txt')
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> # Select a single fold
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> fold_ID = 10
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>> fold_array = -1 * np.ones(npoints)
    >>> fold_array[indices_fold] = 1
    >>>
    >>> fundi, n_fundi, likelihoods = extract_fundi(fold_array, neighbor_lists,
    >>>     depth_file, mean_curv_file, min_curv_vec_file,
    >>>     min_distance=5, thr=0.5, use_only_endpoints=True, compute_local_depth=True)
    >>>
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_extract_fundi.vtk', fundi, 'fundi', folds)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_fundi.vtk')

    """
    import os
    import numpy as np
    from time import time

    from mindboggle.features.fundi import compute_likelihood, find_anchors, connect_points
    from mindboggle.utils.mesh import skeletonize, extract_endpoints
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars

    # Load depth and curvature values from VTK and text files
    faces, lines, indices, points, npoints, depths, \
        name = read_vtk(depth_file, return_first=True, return_array=True)
    points = np.array(points)
    mean_curvatures, name = read_scalars(mean_curvature_file,
                                         return_first=True, return_array=True)
    min_directions = np.loadtxt(min_curvature_vector_file)

    # For each fold region...
    n_folds = len([x for x in np.unique(folds) if x > -1])
    print("Extract a fundus from each of {0} regions...".format(n_folds))
    t1 = time()
    Z = np.zeros(npoints)
    fundi = -1 * np.ones(npoints)
    likelihoods = np.copy(fundi)

    unique_fold_IDs = np.unique(folds)
    unique_fold_IDs = [x for x in unique_fold_IDs if x >= 0]
    count = 0
    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:

            print('  Region {0}:'.format(fold_ID))

            # Compute fundus likelihood values
            local_depths = depths[indices_fold]
            if compute_local_depth:
                local_depths = local_depths / np.max(local_depths)
            fold_likelihoods = compute_likelihood(local_depths,
                                                  mean_curvatures[indices_fold])
            likelihoods[indices_fold] = fold_likelihoods

            # Find fundus points
            fold_indices_anchors = find_anchors(points[indices_fold],
                                                fold_likelihoods,
                                                min_directions[indices_fold],
                                                min_distance, thr)
            indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
            n_anchors = len(indices_anchors)
            if n_anchors > 1:
                t2 = time()
                likelihoods_fold = Z.copy()
                likelihoods_fold[indices_fold] = fold_likelihoods

                # Connect fundus points and extract fundus.
                # If using only endpoints to connect fundus vertices,
                # run in two stages -- fast, to get a rough skeleton
                # to extract the endpoints, then slow, to get fundi
                if not use_only_endpoints:
                    print('    Connect {0} fundus points...'.format(n_anchors))
                    B = connect_points(indices_anchors, indices_fold,
                        likelihoods_fold, neighbor_lists)
                    indices_skeleton = [i for i,x in enumerate(B) if x > 0]
                else:
                    L = likelihoods_fold.copy()
                    L[L > thr] = 1
                    L[L <= thr] = 0
                    S = skeletonize(L, indices_anchors, neighbor_lists)
                    indices_skeleton = [i for i,x in enumerate(S) if x > 0]
                    indices_endpoints = extract_endpoints(indices_skeleton,
                                                          neighbor_lists)
                    indices_endpoints = [x for x in indices_endpoints
                                         if x in indices_anchors]
                    if len(indices_endpoints) > 1:
                        print('    Connect {0} fundus endpoints...'.
                              format(len(indices_endpoints)))
                        B = connect_points(indices_endpoints, indices_fold,
                                           likelihoods_fold, neighbor_lists)
                        indices_skeleton = [i for i,x in enumerate(B) if x > 0]
                    else:
                        indices_skeleton = []

                if len(indices_skeleton) > 1:
                    fundi[indices_skeleton] = fold_ID
                    count += 1
                    print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                          format(n_anchors, time() - t2))

    n_fundi = count
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.format(n_fundi, time() - t1))

    #---------------------------------------------------------------------------
    # Return fundi, number of fundi, likelihood values, and file name
    #---------------------------------------------------------------------------
    if save_file:

        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(depth_file, fundi_file, [fundi, likelihoods],
                        'fundi_likelihoods', fundi)

    else:
        fundi_file = None

    return fundi, n_fundi, likelihoods, fundi_file


# Example
if __name__ == "__main__" :
    import os
    from mindboggle.utils.io_vtk import read_faces_points, read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.features.fundi import extract_fundi
    from mindboggle.utils.mesh import plot_vtk

    path = os.environ['MINDBOGGLE_DATA']

    depth_file = os.path.join(path, 'arno',
                                         'measures', 'lh.pial.depth.vtk')
    labels_file = os.path.join(path, 'arno',
                               'labels', 'lh.labels.DKT25.manual.vtk')
    mean_curvature_file = os.path.join(path, 'arno',
                                         'measures', 'lh.pial.curv.avg.vtk')
    min_curvature_vector_file = os.path.join(path, 'arno',
                                         'measures', 'lh.pial.curv.min.dir.txt')
    sulci_file = os.path.join(path, 'arno',
                                         'features', 'lh.sulci.vtk')

    faces, points, npoints = read_faces_points(depth_file)
    neighbor_lists = find_neighbors(faces, npoints)

    sulci, name = read_scalars(sulci_file, return_first=True, return_array=True)

    fundi, n_fundi, likelihoods = extract_fundi(sulci, neighbor_lists,
        depth_file, mean_curvature_file, min_curvature_vector_file,
        min_distance=5, thr=0.5, use_only_endpoints=True, compute_local_depth=True)

    # Write results to vtk file and view:
    rewrite_scalars(depth_file, 'test_extract_fundi.vtk',
                         [fundi], ['fundi'], sulci)
    plot_vtk('test_extract_fundi.vtk')
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
