#!/usr/bin/python
"""
Connect surface mesh vertices ("anchors").

Connect vertices according to a cost function that penalizes vertices
that do not have high likelihood values and have Hidden Markov Measure Field
(HMMF) values different than their neighbors.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

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


#-----------------------
# Test for simple points
#-----------------------
def simple_test(faces, index, values, thr, neighbors):
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
    neighbors : list of lists of indices of neighboring vertices

    Returns
    -------
    sp : simple point or not?: Boolean
    n_inside : number of neighboring vertices greater than threshold

    """

    import numpy as np

    # Optional:
    # from find_points import find_neighbors

    run_find_neighbors = False
    if run_find_neighbors:
        from find_points import find_neighbors

    # Make sure arguments are numpy arrays
    if type(faces) != np.ndarray:
        faces = np.array(faces)
    if type(values) != np.ndarray:
        values = np.array(values)

    # Find neighbors to the input vertex, and binarize them
    # into those greater than a threshold, thr,
    # and those less than or equal to thr ("inside" and "outside").
    # Count the number of "inside" and "outside" neighbors
    if run_find_neighbors:
        I_neighbors = neighbors
    else:
        I_neighbors = neighbors[index]
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
            if run_find_neighbors:
                new_neighbors = find_neighbors(faces, inside[i_in])
            else:
                new_neighbors = neighbors[inside[i_in]]
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


#===============
# Connect points
#===============
def skeletonize(B, indices_to_keep, faces, neighbor_lists):
    """
    Skeletonize a binary numpy array into 1-vertex-thick curves.

    Parameters
    ----------
    B : binary [#vertices x 1] numpy array
    indices_to_keep : indices to retain
    faces : indices of triangular mesh vertices: [#faces x 3] numpy array
    neighbor_lists : lists of lists of neighboring vertices

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

#===============
# Connect points
#===============
def connect_points(anchors, faces, indices, L, thr, neighbor_lists):
    """
    Connect vertices in a surface mesh to create a curve.

    The goal of this algorithm is to assign each vertex a locally optimal
    Hidden Markov Measure Field (HMMF) value.

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
    neighbor_lists : lists of lists of neighboring vertices
                    (optional: empty list)

    Returns
    -------
    fundus : [#vertices x 1] numpy array

    """

    import numpy as np

    # Optional:
    # from find_points import find_neighbors

    if len(neighbor_lists) == 0:
        from find_points import find_neighbors

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

    # Find neighbors for each vertex
    # (extract_folds() should have found most, if not all, neighbors)
    if len(neighbor_lists) > 0:
        N = neighbor_lists
        #for index in indices:
        #    if not len(N[index]):
        #        N[index] = find_neighbors(faces, index)
    else:
        N = [[] for x in L]
        for index in indices:
            N[index] = find_neighbors(faces, index)

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
                        update, n_in = simple_test(faces, index, H_new, thr, N)
                    elif H[index] <= thr <= H_test:
                        update, n_in = simple_test(faces, index, 1 - H_new, thr, N)

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
            if not np.mod(count, 20):
                print('      {0}: points crossing threshold ({1:.6f}):   {2}'.
                      format(count, delta_cost, delta_points))

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
    skeleton = skeletonize(H, anchors, faces, N)
    print('      Removed {0} points to create one-vertex-thin skeletons'.
          format(n_points - sum(skeleton)))

    return skeleton
