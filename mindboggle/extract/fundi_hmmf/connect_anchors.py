#!/usr/bin/python
"""
Connect surface mesh vertices ("anchors").

Connect vertices according to Hidden Markov Measure Field (HMMF)
probability values.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from find_neighbors import find_neighbors

verbose = 1

#--------------------
# Compute probability
#--------------------
def prob(wt_likelihood, likelihood, wt_neighbors, hmmf, hmmf_neighbors):
    """
    Compute Hidden Markov Measure Field (HMMF) probability for a given vertex.

    p = hmmf * np.sqrt((wt_likelihood - likelihood)**2)
        + wt_neighbors * sum((hmmf - hmmf_neighbors)**2)

    term 1 promotes high likelihood values
    term 2 promotes smoothness of the HMMF values

    Inputs:
    ------
    wt_likelihood: influence of likelihood on probability (term 1): float > 1
    wt_neighbors: weight influence of neighbors on probability (term 2)
    likelihood: likelihood value in interval [0,1]
    hmmf: HMMF value
    hmmf_neighbors: HMMF values of neighboring vertices

    Output:
    ------
    probability

    """
    return hmmf * (wt_likelihood - likelihood) + \
           wt_neighbors * sum((hmmf - hmmf_neighbors)**2)

#-----------------------
# Test for simple points
#-----------------------
def simple_test(faces, index, values, thr, neighbors):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when removed from or added to an object
    (e.g., a curve) on a surface mesh does not alter the topology of the object.

    Inputs:
    ------
    faces: [#faces x 3] numpy array
    index: index of vertex
    values: values: [#vertices x 1] numpy array
    thr: threshold
    neighbors: list of lists of indices of neighboring vertices

    Output:
    ------
    sp: simple point or not? (1 or 0)

    Calls:
    -----
    find_neighbors()  (optional)

    """

    # Find neighbors to the input vertex, and binarize them
    # into those greater than the likelihood threshold, thr,
    # and those less than or equal to thr ("inside" and "outside").
    # Count the number of "inside" and "outside" neighbors
    nlist = 1
    if nlist:
        I_neighbors = neighbors[index]
    else:
        I_neighbors = neighbors
    neighbor_values = values[I_neighbors]
    inside = [I_neighbors[i] for i,x in enumerate(neighbor_values) if x > thr]
    n_inside = len(inside)
    n_outside = len(I_neighbors) - n_inside

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point
    if n_outside * n_inside == 0:
        sp = 0
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point
    elif n_outside == 1 or n_inside == 1:
        sp = 1
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point
    else:
        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' indices in a sublist of "N"
        labels = range(1, n_inside + 1)
        N = []
        for i_in in range(n_inside):
            if nlist:
                new_neighbors = neighbors[inside[i_in]]
            else:
                new_neighbors = find_neighbors(faces, inside[i_in])
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
            sp = 1
        else:
            sp = 0

    return sp

#========================
# Connect anchor vertices
#========================
def connect_anchors(anchors, faces, indices, L, thr, neighbor_lists):
    """
    Connect anchor vertices in a surface mesh to create a curve.

    ????[Explain HMMF approach]

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
    Parameters for computing the probabilities and probability gradients:
      wt_likelihood: weight influence of likelihood on probability
      wt_neighbors: weight influence of neighbors on probability
      decrement: the amount that the HMMF values are decremented
    Parameters to speed up optimization and terminate the algorithm:
      min_C: minimum value to fix probabilities at very low values
      min_change: minimum change in the sum of probabilities
      n_tries_no_change: #times the loop can continue even without any change
      max_count: maximum #iterations

    Output:
    ------
    fundus: [#vertices x 1] numpy array

    Calls:
    -----
    find_neighbors()
    prob()
    simple_test()

    """

    #-----------
    # Parameters
    #-------------------------------------------------------------------------
    # prob() and probability gradient parameters
    wt_likelihood = 1.1  # float > 1: affects likelihood's influence on probability
    wt_neighbors = 0.4  # weight influence of neighbors on probability
    decrement = 0.05  # the amount that the HMMF values are decremented

    # Parameters to speed up optimization and for termination of the algorithm
    min_C = 0.01  # minimum value to fix probabilities at low values
    min_change = 0.0001  # minimum change in the sum of probabilities
    n_tries_no_change = 3  # #times loop can continue even without any change
    max_count = 100  # maximum number of iterations
    #-------------------------------------------------------------------------

    # Initialize all likelihood values within sulcus as greater than 0.5
    # and less than or equal to 1.0.
    # This is necessary to guarantee correct topology
    L_init = L.copy() + thr + 0.000001
    L_init[L == 0.0] = 0
    L_init[L_init > 1.0] = 1

    # Initialize array of connected vertices C
    # with likelihood values greater than thr and a 1 for each anchor vertex
    n_vertices = len(indices)
    Z = np.zeros(len(L))
    C = Z.copy()
    C[L_init > thr] = L_init[L_init > thr]
    C[anchors] = 1

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

    # Assign probability values to each vertex
    probs = Z.copy()
    probs[indices] = [prob(wt_likelihood, L[i], wt_neighbors, C[i], C[N[i]])
                      for i in indices]
    print('      Assigned probability values to {} vertices'.
          format(n_vertices))

    # Loop until count reaches max_count or until end_flag equals zero
    # (end_flag is used to allow the loop to continue even if there is
    #  no change for n_tries_no_change times)
    count = 0
    end_flag = 0
    Cnew = C.copy()
    while end_flag < n_tries_no_change and count < max_count:

        # Assign each vertex a locally optimal HMMF value (C[i])
        for i in indices:

          # Do not update anchor point probabilities
          if i not in anchors:

            # Continue if HMMF value is greater than C_threshold
            # (to fix when at very low values, to speed up optimization)
            if C[i] > min_C:

                # Compute the probability gradient for the HMMF value
                q = max([C[i] - decrement, 0])
                prob_decr = prob(wt_likelihood, L[i],
                                 wt_neighbors, q, C[N[i]])
                test_value = C[i] - (probs[i] - prob_decr)

                # Update the HMMF value if just above the threshold
                # such that the decrement makes it cross the threshold
                # and the vertex is a "simple point"
                # (Cnew[i] is not changed yet since simple_test()
                #  only considers its neighbors)
                if C[i] >= thr >= test_value:
                    update = simple_test(faces, i, Cnew, thr, N)
                elif C[i] <= thr <= test_value:
                    update = simple_test(faces, i, 1 - Cnew, thr, N)
                # Update the HMMF value if far from the threshold
                else:
                    update = 1

                # Update the HMMF and probability values
                if update:
                    if test_value < 0:
                        test_value = 0.0
                    elif test_value > 1:
                        test_value = 1.0
                    Cnew[i] = test_value
                    probs[i] = prob(wt_likelihood, L[i],
                                    wt_neighbors, Cnew[i], C[N[i]])

        # Sum the probability values across all vertices
        # and tally the number of HMMF values with probability > thr.
        # After iteration 1, compare current and previous values.
        # If the values are similar, increment end_flag.
        sum_probs = sum(probs)
        n_points = sum([1 for x in C if x > thr])
        if verbose == 2:
            print('        {} vertices...'.format(n_points))

        if count > 0:
            if n_points == n_points_previous:
                if (sum_probs_previous - sum_probs) / n_vertices < min_change:
                    end_flag += 1

        # Reset for next iteration
        sum_probs_previous = sum_probs
        n_points_previous = n_points
        C = Cnew

        count += 1

    print('      Updated HMMF values')

    return C.tolist()
