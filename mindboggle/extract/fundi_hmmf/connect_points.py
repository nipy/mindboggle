#!/usr/bin/python
"""
Connect surface mesh vertices ("anchors").

Connect vertices according to a cost function that penalizes vertices
that do not have high likelihood values and have Hidden Markov Measure Field
(HMMF) values different than their neighbors.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from compute_values import compute_cost

# Optional in simple_test() and in connect_points():
# from find_points import find_neighbors


#-----------------------
# Test for simple points
#-----------------------
def simple_test(faces, index, values, thr, neighbors):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when added to or removed from an object
    (e.g., a curve) on a surface mesh does not alter the object's topology.

    Inputs:
    ------
    faces: [#faces x 3] numpy array
    index: index of vertex
    values: values: [#vertices x 1] numpy array
    thr: threshold
    neighbors: list of lists of indices of neighboring vertices

    Output:
    ------
    sp: simple point or not?: Boolean

    Calls:
    -----
    find_neighbors()  [optional]

    """

    run_find_neighbors = False
    if run_find_neighbors:
        from find_points import find_neighbors

    # Find neighbors to the input vertex, and binarize them
    # into those greater than the likelihood threshold, thr,
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

    return sp


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
    threshold such that a step_down makes it cross the threshold,
    and the vertex is a "simple point" (its addition/removal alters topology).

    Inputs:
    ------
    anchors: list of indices of vertices to connect (should contain >=2)
    faces: indices of triangular mesh vertices: [#faces x 3] numpy array
    indices: list of indices of vertices
    L: likelihood values: [#vertices in mesh x 1] numpy array
    thr: likelihood threshold
    neighbor_lists: lists of lists of neighboring vertices (optional: empty list)

    Parameters:
    ----------
    Parameters for computing the cost and cost gradients:
      wL: weight influence of likelihood on the cost function
      wN: weight influence of neighbors on the cost function
      step_down: the amount that the HMMF values are step_downed
    Parameters to speed up optimization and terminate the algorithm:
      min_H: minimum HMMF value to fix very low values
      min_change: minimum change in the sum of costs
      n_tries_no_change: #times the loop can continue even without any change
      max_count: maximum #iterations

    Output:
    ------
    fundus: [#vertices x 1] numpy array

    Calls:
    -----
    find_neighbors() [optional]
    compute_cost()
    simple_test()

    """

    if len(neighbor_lists) == 0:
        from find_points import find_neighbors

    #-----------
    # Parameters
    #-------------------------------------------------------------------------
    # Cost and cost gradient parameters
    wL = 1 #1.1 # 1.1  # weight of likelihood on cost function
    wN = 0.5 # 0.4 # 0.4  # initial weight of neighbors on cost function
    step_down = 0.05 # the amount that HMMF values are stepped down

    # Parameters to speed up optimization and for termination of the algorithm
    gradient_init = 0.1  # initialize gradient factor
    gradient_step = 0.001  # step gradient factor up each iteration
    min_cost_change = 0.000001  # minimum change in the sum of costs
    n_tries_no_change = 3  # sequential loops without sufficient change
    max_count = 1000  # maximum number of iterations (in case no convergence)
    #-------------------------------------------------------------------------

    # Initialize all Hidden Markov Measure Field (HMMF) values with
    # likelihood values (except 0) normalized to the interval (0.5, 1.0]
    # (to guarantee correct topology) and take those values that are greater
    # than the likelihood threshold.  Assign a 1 for each anchor point.
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
    C[indices] = [compute_cost(wL, wN, L[i], H[i], H[N[i]])[0] for i in indices]

    # Loop until count reaches max_count or until end_flag equals zero
    # (end_flag is used to allow the loop to continue even if there is
    #  no change for n_tries_no_change times)
    count = 0
    end_flag = 0
    H_new = H.copy()
    while end_flag < n_tries_no_change and count < max_count:

        # At each iteration, step up gradient factor
        gradient_factor = gradient_init + gradient_step * count

        # At each iteration, de-weight neighborhood influence of cost function
#        wN = wN_init - wN_step * count

        # For each index
        for index in indices:

            # Do not update anchor point costs
            if index not in anchors:

                # Continue if the HMMF value is greater than the gradient step
                # (to fix when at very low values, to speed up optimization)
                if H[index] > gradient_step:

                    # Compute the cost gradient for the HMMF value
                    H_down = max([H[index] - step_down, 0])

                    cost_down,cw,cn = compute_cost(wL, wN, L[index], H_down, H[N[index]])

                    #print(cost_down, cw, cn)

                    H_test = H[index] - gradient_factor * (C[index] - cost_down)

                    # Update the HMMF value if near the threshold
                    # such that a step_down makes it cross the threshold,
                    # and the vertex is a "simple point"
                    # Note: H_new[index] is not changed yet since simple_test()
                    #       only considers its neighbors
                    if H[index] >= thr >= H_test:
                        update = simple_test(faces, index, H_new, thr, N)
                    elif H[index] <= thr <= H_test:
                        update = simple_test(faces, index, 1 - H_new, thr, N)

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
                        C[index] = compute_cost(wL, wN, L[index],
                                                H_new[index], H[N[index]])[0]

        # Sum the cost values across all vertices and tally the number
        # of HMMF values greater than the threshold.
        # After iteration 1, compare current and previous values.
        # If the values are similar, increment end_flag.
        sum_C = sum(C)
        n_points = sum([1 for x in H if x > thr])

        # Terminate the loop if there are insufficient changes
        if count > 0:
            print('      {}: factor={:.3f}; -points={}; delta cost={:.8f}'.
                  format(count, gradient_factor, n_points_previous - n_points,
                         (sum_C_previous - sum_C) / n_vertices))
            if n_points == n_points_previous:
                if (sum_C_previous - sum_C) / n_vertices < min_cost_change:
                    end_flag += 1
            else:
                end_flag = 0

        # Reset for next iteration
        sum_C_previous = sum_C
        n_points_previous = n_points
        H = H_new

        count += 1

    print('      Updated hidden Markov measure field (HMMF) values')

    H_binary = H.copy()

    return H.tolist(), H_binary.tolist()
