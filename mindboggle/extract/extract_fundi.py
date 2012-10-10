#!/usr/bin/python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import numpy as np
#from time import time
#from utils.mesh_operations import simple_test, skeletonize, find_anchors,
#                                  extract_endpoints
#from mindboggle.extract.extract_fundi import compute_likelihood, connect_points
#from mindboggle.utils.mesh_operations import find_anchors, extract_endpoints,
#                                             simple_test, skeletonize
#from mindboggle.utils.io_vtk import load_scalar

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
    depths : depth values in [0,1]: [#sulcus vertices x 1] numpy array
    curvatures : mean curvature values in [-1,1] [#sulcus vertices x 1] array

    Returns
    -------
    likelihoods : likelihood values [#sulcus vertices x 1] numpy array

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if type(depths) != np.ndarray:
        depths = np.array(depths)
    if type(curvatures) != np.ndarray:
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
    hmmf_neighbors : numpy array
        HMMF values of neighboring vertices
    wN : float
        weight influence of neighbors on cost (term 2)

    Returns
    -------
    cost : ``float``

    """

#    # Make sure argument is a numpy array
#    import numpy as np
#    if type(hmmf_neighbors) != np.ndarray:
#        hmmf_neighbors = np.array(hmmf_neighbors)

    if len(hmmf_neighbors):
        cost = hmmf * (1 - likelihood) +\
               wN * sum(abs(hmmf - hmmf_neighbors)) / len(hmmf_neighbors)
    else:
        exit('ERROR: No HMMF neighbors to compute cost.')

    return cost

#---------------
# Connect points
#---------------
def connect_points(anchors, faces, indices, L, neighbor_lists):
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

        ``H_step``: the amount that the HMMF values are H_steped

    Parameters to speed up optimization and terminate the algorithm:

        ``min_H``: minimum HMMF value to fix very low values

        ``min_change``: minimum change in the sum of costs

        ``n_tries_no_change``: #times the loop can continue even without any change

        ``max_count``: maximum #iterations

    Parameters
    ----------
    anchors : list of indices of vertices to connect (should contain > 1)
    faces : indices of triangular mesh vertices: [#faces x 3] numpy array
    indices : list of indices of vertices
    L : likelihood values: [#vertices in mesh x 1] numpy array
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    skeleton : [#vertices x 1] numpy array

    Example
    -------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors, find_anchors
    >>> from mindboggle.extract.extract_fundi import connect_points, compute_likelihood
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.avg.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.min.dir.txt')
    >>> sulci_file = os.path.join(data_path, 'results', 'features',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, True)
    >>> points, faces, mean_curvatures, n_vertices = load_scalar(mean_curvature_file,
    >>>                                                          True)
    >>> points, faces, sulcus_IDs, n_vertices = load_scalar(sulci_file, True)
    >>> min_directions = np.loadtxt(min_curvature_vector_file)
    >>> sulcus_ID = 1
    >>> indices_fold = [i for i,x in enumerate(sulcus_IDs) if x == sulcus_ID]
    >>> fold_likelihoods = compute_likelihood(depths[indices_fold],
    >>>                                       mean_curvatures[indices_fold])
    >>> likelihoods = np.zeros(len(points))
    >>> likelihoods[indices_fold] = fold_likelihoods
    >>> fold_indices_anchors = find_anchors(points[indices_fold, :],
    >>>     fold_likelihoods, min_directions[indices_fold], 5, 0.5)
    >>> indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
    >>> likelihoods_fold = np.zeros(len(points))
    >>> likelihoods_fold[indices_fold] = fold_likelihoods
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> H = connect_points(indices_anchors, faces, indices_fold,
    >>>     likelihoods_fold, neighbor_lists)
    >>> # Write results to vtk file and view with mayavi2:
    >>> H[indices_anchors] = 2
    >>> rewrite_scalars(depth_file, 'test_connect_points.vtk',
    >>>                 H, H)
    >>> os.system('mayavi2 -m Surface -d test_connect_points.vtk &')

    """
    import numpy as np
    from mindboggle.utils.mesh_operations import simple_test, skeletonize

    # Make sure arguments are numpy arrays
    if type(faces) != np.ndarray:
        faces = np.array(faces)
    if type(L) != np.ndarray:
        L = np.array(L)

    #-----------
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
    n_vertices = len(indices)
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
    print(len(indices))
    print(np.unique(H))
    print(np.unique(L))
    print(L[0], H[0], H[N[0]], wN_max)
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
                # (to fix when at very low values, to speed up optimization)
                if H[index] > H_min:

                    # Compute the cost gradient for the HMMF value
                    H_down = max([H[index] - H_step, 0])
                    cost_down = compute_cost(L[index], H_down, H[N[index]], wN)
                    H_test = H[index] - gradient_factor * (C[index] - cost_down)

                    # Update the HMMF value if near the threshold
                    # such that a step makes it cross the threshold,
                    # and the vertex is a "simple point"
                    # Note: H_new[index] is not changed yet since
                    #       simple_test() only considers its neighbors
                    if H[index] >= 0.5 >= H_test:
                        update, n_in = simple_test(index, H_new, N)
                    elif H[index] <= 0.5 <= H_test:
                        update, n_in = simple_test(index, 1 - H_new, N)

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
            delta_cost = (costs_previous - costs) / n_vertices
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
def extract_fundi(fold_IDs, neighbor_lists, depth_file,
                  mean_curvature_file, min_curvature_vector_file,
                  min_distance=5, thr=0.5, use_only_endpoints=True):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Parameters
    ----------
    fold_IDs : list or numpy array
        fold IDs (default = -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    depth_file : str
        surface mesh file in VTK format with faces and scalar values
    mean_curvature_file : str
        surface mesh file in VTK format with scalar values
    min_curvature_vector_file : str
        surface mesh file in VTK format with scalar values
    min_fold_size : int
        minimum fold size (number of vertices)
    min_distance :  int
        minimum distance
    thr :  float
        likelihood threshold
    use_only_endpoints : Boolean
        use endpoints to construct fundi (or all anchor points)?

    Returns
    -------
    fundus_IDs : array of integers
        fundus IDs for all vertices, with -1s for non-fundus vertices
    n_fundi :  int
        number of sulcus fundi

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.extract.extract_fundi import extract_fundi
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.avg.vtk')
    >>> min_curvature_vector_file = os.path.join(data_path, 'measures',
    >>>     '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.min.dir.txt')
    >>> sulci_file = os.path.join(data_path, 'results', 'features',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, True)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> points, faces, sulcus_IDs, n_vertices = load_scalar(sulci_file, True)
    >>> fundus_IDs, n_fundi = extract_fundi(sulcus_IDs, neighbor_lists,
    >>>     depth_file, mean_curvature_file, min_curvature_vector_file,
    >>>     min_distance=5, thr=0.5, use_only_endpoints=True)
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalars(depth_file, 'test_extract_fundi.vtk',
    >>>                 fundus_IDs, fundus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_fundi.vtk &')
    >>> # Write and view manual labels restricted to sulci:
    >>> rewrite_scalars(depth_file, 'test_extract_sulci_labels.vtk',
    >>>                 labels, sulcus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci_labels.vtk &')

    """
    import numpy as np
    from time import time

    from mindboggle.extract.extract_fundi import compute_likelihood, connect_points
    from mindboggle.utils.mesh_operations import find_anchors, skeletonize, extract_endpoints
    from mindboggle.utils.io_vtk import load_scalar

    # Load depth and curvature values from VTK and text files
    points, faces, depths, n_vertices = load_scalar(depth_file, True)
    points, faces, mean_curvatures, n_vertices = load_scalar(mean_curvature_file,
                                                             True)
    min_directions = np.loadtxt(min_curvature_vector_file)

    # For each fold region...
    n_folds = len([x for x in np.unique(fold_IDs) if x > -1])
    print("Extract a fundus from each of {0} regions...".format(n_folds))
    t1 = time()
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()
    fundus_IDs = -1 * np.ones(n_vertices)

    unique_fold_IDs = np.unique(fold_IDs)
    unique_fold_IDs = [x for x in unique_fold_IDs if x >= 0]
    count = 0
    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(fold_IDs) if x == fold_ID]
        if len(indices_fold):

            print('  Region {0}:'.format(fold_ID))

            # Compute fundus likelihood values
            fold_likelihoods = compute_likelihood(depths[indices_fold],
                                                  mean_curvatures[indices_fold])
            likelihoods[indices_fold] = fold_likelihoods

            # Find fundus points
            fold_indices_anchors = find_anchors(points[indices_fold, :],
                                                fold_likelihoods,
                                                min_directions[indices_fold],
                                                min_distance, thr)
            indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
            n_anchors = len(indices_anchors)
            if n_anchors > 1:

                # Connect fundus points and extract fundus
                print('    Connect {0} fundus points...'.format(n_anchors))
                t2 = time()
                likelihoods_fold = Z.copy()
                likelihoods_fold[indices_fold] = fold_likelihoods

                # If using only endpoints to connect fundus vertices,
                # run in two stages -- fast, to get a rough skeleton
                # to extract the endpoints, then slow, to get fundi
                if not use_only_endpoints:
                    B = connect_points(indices_anchors, faces, indices_fold,
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
                        B = connect_points(indices_endpoints, faces, indices_fold,
                                           likelihoods_fold, neighbor_lists)
                        indices_skeleton = [i for i,x in enumerate(B) if x > 0]
                    else:
                        indices_skeleton = []

                if len(indices_skeleton) > 1:
                    fundus_IDs[indices_skeleton] = fold_ID
                    count += 1
                    print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                          format(n_anchors, time() - t2))

    n_fundi = count
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.format(n_fundi, time() - t1))

    return fundus_IDs, n_fundi


# Example
if __name__ == "__main__" :
    import os
    from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    from mindboggle.utils.mesh_operations import find_neighbors
    from mindboggle.extract.extract_fundi import extract_fundi

    data_path = os.environ['MINDBOGGLE_DATA']

    depth_file = os.path.join(data_path, 'measures',
                 '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')

    mean_curvature_file = os.path.join(data_path, 'measures',
        '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.avg.vtk')

    min_curvature_vector_file = os.path.join(data_path, 'measures',
        '_hemi_lh_subject_MMRR-21-1', 'lh.pial.curv.min.dir.txt')

    sulci_file = os.path.join(data_path, 'results', 'features',
                 '_hemi_lh_subject_MMRR-21-1', 'sulci.vtk')

    points, faces, depths, n_vertices = load_scalar(depth_file, True)
    neighbor_lists = find_neighbors(faces, len(points))

    points, faces, sulcus_IDs, n_vertices = load_scalar(sulci_file, True)

    fundus_IDs, n_fundi = extract_fundi(sulcus_IDs, neighbor_lists,
        depth_file, mean_curvature_file, min_curvature_vector_file,
        min_distance=5, thr=0.5, use_only_endpoints=True)

    # Write results to vtk file and view with mayavi2:
    rewrite_scalars(depth_file, 'test_extract_fundi.vtk',
                    fundus_IDs, fundus_IDs)
    os.system('mayavi2 -m Surface -d test_extract_fundi.vtk &')
