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
#from utils.mesh_operations import simple_test, skeletonize

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
    likelihood : likelihood value in interval [0,1]
    hmmf : HMMF value
    hmmf_neighbors : HMMF values of neighboring vertices: numpy array
    wN : weight influence of neighbors on cost (term 2)

    Returns
    -------
    cost : ``float``

    """

#    # Make sure argument is a numpy array
#    import numpy as np
#    if type(hmmf_neighbors) != np.ndarray:
#        hmmf_neighbors = np.array(hmmf_neighbors)

    cost = hmmf * (1 - likelihood) +\
           wN * sum(abs(hmmf - hmmf_neighbors)) / len(hmmf_neighbors)

    return cost

#---------------
# Connect points
#---------------
def connect_points(anchors, faces, indices, L, thr, neighbor_lists):
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
    anchors : list of indices of vertices to connect (should contain >=2)
    faces : indices of triangular mesh vertices: [#faces x 3] numpy array
    indices : list of indices of vertices
    L : likelihood values: [#vertices in mesh x 1] numpy array
    thr : likelihood threshold
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    fundus : [#vertices x 1] numpy array

    """
    import numpy as np
    from utils.mesh_operations import simple_test, skeletonize

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
    H_min = thr - H_step  # minimum HMMF value to be processed

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
    n_vertices = len(indices)
    C = np.zeros(len(L))
    H = C.copy()
    H_init = (L + 1.000001) / 2
    H_init[L == 0.0] = 0
    H_init[H_init > 1.0] = 1
    H[H_init > thr] = H_init[H_init > thr]
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
                    if H[index] >= thr >= H_test:
                        update, n_in = simple_test(index, H_new, thr, N)
                    elif H[index] <= thr <= H_test:
                        update, n_in = simple_test(index, 1 - H_new, thr, N)

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
        n_points = sum([1 for x in H if x > thr])

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
                print('      Iteration {0}: {2} points crossing threshold'.
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
    H[H > thr] = 1
    H[H <= thr] = 0
    n_points = sum(H)

    # Skeletonize
    skeleton = skeletonize(H, anchors, N)
    print('      Removed {0} points to create one-vertex-thin skeletons'.
          format(int(n_points - sum(skeleton))))

    return skeleton

#==================
# Extract all fundi
#==================
def extract_fundi(folds, n_folds, neighbor_lists,
                  depth_file, mean_curvature_file, min_curvature_vector_file,
                  min_fold_size=50, min_distance=5, thr=0.5):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Parameters
    ----------
    folds : list or numpy array
        fold IDs (default = -1)
    n_folds :  int
        number of folds
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

    Returns
    -------
    fundi :  numpy array of fundi

    """
    import numpy as np
    from time import time

    from extract.extract_fundi import compute_likelihood, connect_points
    from utils.mesh_operations import find_anchors
    from utils.io_vtk import load_scalar

    # Convert folds array to a list of lists of vertex indices
    index_lists_folds = [np.where(folds == i)[0].tolist()
                         for i in range(n_folds)]

    # Load depth and curvature values from VTK and text files
    vertices, Faces, depths, n_vertices = load_scalar(depth_file, return_arrays=1)
    vertices, Faces, mean_curvatures, n_vertices = load_scalar(mean_curvature_file,
                                                               return_arrays=1)
    min_directions = np.loadtxt(min_curvature_vector_file)

    # For each fold...
    print("Extract a fundus from each of {0} folds...".format(n_folds))
    t1 = time()
    fundus_lists = []
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()

    for i_fold, indices_fold in enumerate(index_lists_folds):

        print('  Fold {0} of {1}:'.format(i_fold + 1, n_folds))

        # Compute fundus likelihood values
        fold_likelihoods = compute_likelihood(depths[indices_fold],
                                              mean_curvatures[indices_fold])
        likelihoods[indices_fold] = fold_likelihoods

        # If the fold has enough high-likelihood vertices, continue
        likelihoods_thr = sum(fold_likelihoods > thr)
        print('    {0} vertices with fundus likelihood value > {1}'.
              format(likelihoods_thr, thr)) #min_fold_size
        if likelihoods_thr > min_fold_size:

            # Find fundus points
            fold_indices_anchors = find_anchors(vertices[indices_fold, :],
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

                H = connect_points(indices_anchors, Faces, indices_fold,
                                   likelihoods_fold, thr, neighbor_lists)
                fundus_lists.append(H.tolist())
                print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                      format(n_anchors, time() - t2))
            else:
                fundus_lists.append([])
        else:
            fundus_lists.append([])

    fundi = -1 * np.ones(n_vertices)
    count = 0
    for fundus in fundus_lists:
        if len(fundus) > 0:
            fundi[fundus] = count
            count += 1

    print('  ...Extracted fundi ({0:.2f} seconds)'.format(time() - t1))

    return fundi
