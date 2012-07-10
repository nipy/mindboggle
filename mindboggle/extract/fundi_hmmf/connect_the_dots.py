#!/usr/bin/python
"""
Connect surface mesh vertices ("dots") according to Hidden Markov Measure Field
(HMMF) probability values.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from find_neighbors import find_neighbors

#--------------------
# Compute probability
#--------------------
def prob(wt_likelihood, likelihood, wt_neighbors, hmmf, hmmf_neighbors):
    """
    Compute Hidden Markov Measure Field (HMMF) probability for a given vertex.

    p = -hmmf * np.sqrt(wt_likelihood - likelihood)**2
        - wt_neighbors * sum(hmmf - hmmf_neighbors)**2

    term 1 promotes high likelihood values
    term 2 promotes smoothness of the HMMF values

    Inputs:
    ------
    wt_likelihood: weight influence of likelihood on probability (term 1)
    wt_neighbors: weight influence of neighbors on probability (term 2)
    likelihood: likelihood value
    hmmf: HMMF value
    hmmf_neighbors: HMMF values of neighboring vertices

    Output:
    ------
    p: probability

    """
    p = -hmmf * np.sqrt(wt_likelihood - likelihood)**2 - \
        wt_neighbors * sum(hmmf - hmmf_neighbors)**2
    return p

#-----------------------
# Test for simple points
#-----------------------
def simple_point_test(faces, index, values, thr, max_neighbors=20):
    """
    Test to see if vertex is a "simple point".

    A simple point is a vertex that when removed from or added to an object
    (such as a curve) on a surface mesh does not alter the topology of the object.

    Inputs:
    ------
    faces: [#faces x 3] numpy array
    index: index of vertex
    values: likelihood values: [#vertices x 1] numpy array
    thr: likelihood threshold
    max_neighbors: maximum number of neighbors considered for a given vertex

    Output:
    ------
    sp: simple point or not? (1 or 0)

    Calls:
    -----
    find_neighbors()

    """

    # Find neighbors to the input vertex, and binarize them
    # into those greater than the likelihood threshold, thr,
    # and those less than or equal to thr ("inside" and "outside").
    # Count the number of "inside" and "outside" neighbors
    neighs = find_neighbors(faces, index)
    neighvals = values[neighs]
    neighvals_thr = neighvals > thr
    n_sets = sum(neighvals_thr)
    n_outside = sum(neighvals <= thr)

    # If the number of inside or outside neighbors is zero,
    # than the vertex IS NOT a simple point
    if n_outside * n_sets == 0:
        sp = 0
    # Or if either the number of inside or outside neighbors is one,
    # than the vertex IS a simple point
    elif n_outside == 1 or n_sets == 1:
        sp = 1
    # Otherwise, test to see if all of the inside neighbors share neighbors
    # with each other, in which case the vertex IS a simple point
    else:
        #neighvals_list = [n for n in neighvals > thr]
        neighs_inside = neighs[neighvals_thr]

        # Initialize "sets" numpy array [#neighbors>thr x max_neighbors]
        sets = np.zeros((n_sets, max_neighbors))
        set_labels = range(1, n_sets + 1)

        # Reset vertex value????
        values[index] = -1

        # For each neighbor exceeding the threshold,
        # find its neighbors that also exceed the threshold,
        # then store these neighbors' values in a row of the "sets" array
        for i in range(n_sets):
            current_neighs = find_neighbors(faces, neighs_inside[i])
            current_neighs = current_neighs[values[current_neighs] > thr]
            len_current = len(current_neighs)
            if len_current > max_neighbors - 1:
                len_current = max_neighbors - 1
                current_neighs = current_neighs[range(len_current)]
            sets[i, range(len_current)] = current_neighs
            sets[i, len_current] = neighs_inside[i]

        # Consolidate labels of connected vertices:
        # Loop through neighbors (rows of the "sets" array),
        # reassigning the labels for the rows until each label's
        # row(s) has a unique set of vertices
        change = 1
        while change > 0:
            change = 0

            # Loop through pairs of rows of the "set" array
            # and continue if their two labels are different
            for i in range(n_sets - 1):
                for j in range(i + 1, n_sets):
                    if set_labels[i] != set_labels[j]:

                        # See if the two rows share a (non-zero) vertex
                        unique_i = np.unique(sets[i, :])
                        unique_j = np.unique(sets[j, :])
                        unique_i = unique_i[unique_i > 0]
                        unique_j = unique_j[unique_j > 0]
                        len_union = len(np.union1d(unique_i, unique_j))

                        # Assign the two rows the same label
                        # if they share at least one vertex,
                        # and continue looping
                        if len_union < sum(unique_i > 0) + sum(unique_j > 0):
                            set_labels[i] = set_labels[j]
                            change = 1

        # The vertex is a simple point if all of its neighbors
        # (if any) share neighbors with each other
        n_separate_sets = len(np.unique(set_labels))
        if n_separate_sets == 1:
            sp = 1
        else:
            sp = 0

    return sp

#=================
# Connect the dots
#=================
def connect_the_dots(L, faces, dots, neighbors, indices, thr):
    """
    Connect dots (sulcus vertices in a surface mesh) to create a (fundus) curve.

    ????[Explain HMMF approach]

    Inputs:
    ------
    L: likelihood values for a subset of vertices (zero for the rest):
       [#vertices x 1] numpy array
    faces: [#faces x 3] numpy array
    dots: [#vertices x 1] numpy array
    neighbors: indices of neighbors to each vertex in vertices_subset:
               list of #vertices_subset lists
    indices: indices of subset of vertices: [#vertices_subset x 1] numpy array
    thr: likelihood threshold

    Parameters:
    ----------
    prob() parameters:
      wt_likelihood: weight influence of likelihood on probability
      wt_neighbors: weight influence of neighbors on probability
    Probability gradient parameters:
      step_size
      mult_init
    Loop parameters:
      n_tries_no_change: number of times the loop can continue even without any change
      max_count: maximum number of iterations
    Parameters to speed up optimization:
      diff_threshold
      mult_incr
      C_threshold: a minimum HMMF value to fix probabilities at very low values

    Output:
    ------
    fundus: [#vertices x 1] numpy array

    Calls:
    -----
    prob()
    simple_point_test()

    """

    #-----------
    # Parameters
    #-------------------------------------------------------------------------------------
    # prob() parameters
    wt_likelihood = 1.1  # weight influence of likelihood on probability
    wt_neighbors = 0.4  # weight influence of neighbors on probability

    # Probability gradient parameters
    step_size = 0.05
    mult_init = 0.02

    # Loop parameters
    n_tries_no_change = 3  # number of times the loop can continue even without any change
    max_count = 350  # maximum number of iterations

    # Parameters to speed up optimization
    diff_threshold = 0.0001
    mult_incr = 0.001
    C_threshold = 0.01  # a minimum HMMF value to fix probabilities at very low values
    #-------------------------------------------------------------------------------------

    # Initialize all likelihood values within sulcus greater than 0.5
    # and less than or equal to 1.0.
    # This is necessary to guarantee correct topology
    L_init = (L + 1.000001) / 2.0
    L_init[L == 0] = 0
    L_init[L_init > 1] = 1

    # Initialize array of connected dots,
    # containing likelihood values greater than thr,
    # and 1 for dots greater than thr
    n_vertices = len(L)
    C = np.zeros(n_vertices)
    C[L_init > thr] = L_init[L_init > thr]
    C[dots > thr] = 1

    # Continue if there are at least two candidate dots
    if sum(C > 0) >= 2:
        print(str(sum(C > 0)) + ' initial candidate vertices')

        # Initialize new arrays of connected points and probability values
        Cnew = C.copy()
        probs = np.zeros(n_vertices)
        n_indices = len(indices)

        # Assign probability values to each neighboring vertex
        for i in indices:
            probs[i] = prob(wt_likelihood, L[i], wt_neighbors, C[i], C[neighbors[i]])

        # Loop until count reaches max_count or until end_flag equals zero
        # (end_flag is used to allow the loop to continue even if there is
        #  no change for n_tries_no_change times)
        count = 0
        end_flag = 0
        while end_flag < n_tries_no_change and count < max_count:

            # Define a factor to multiply the probability gradient that will
            # increase the time-step size toward the end of the optimization
            mult = mult_init + count * mult_incr

            # Loop through vertices with positive likelihood values
            # and assign each a locally optimal HMMF value (C[i]) 
            for i in indices:
                # Continue if HMMF value is greater than C_threshold
                # (to fix at very low values, to make the optimization faster)
                if C[i] > C_threshold:

                    # Compute the (negative) probability gradient for the HMMF value
                    q = max([C[i] - step_size, 0])
                    prob_decr = prob(wt_likelihood, L[i], wt_neighbors, q, C[neighbors[i]])
                    decr = mult * (prob_decr - probs[i]) / step_size

                    # Test to update the HMMF value for positive and negative decrements:
                    if decr > 0:
                        # Update the HMMF value if it is close enough to the threshold
                        # that the decrement makes it cross the threshold,
                        # the dot value is below the threshold,
                        # and the vertex is a "simple point"
                        if C[i] > thr >= C[i] - decr:
                            if dots[i] > thr:
                                update = 0
                            else:
                                Cnew_copy = Cnew.copy()
                                Cnew_copy[i] = C[i] - decr
                                update = simple_point_test(faces, i, Cnew_copy, thr)
                        # Or update the HMMF value if it is far from the threshold
                        else:
                            update = 1
                        # Update the HMMF and probability values
                        if update:
                            Cnew[i] = max([C[i] - decr, 0])
                            probs[i] = prob(wt_likelihood, L[i], wt_neighbors,
                                            Cnew[i], C[neighbors[i]])
                    else:
                        # Or if the decrement is negative,
                        # update the HMMF value if it is close enough to the threshold
                        # that the decrement makes it cross the threshold,
                        # and the vertex is a "simple point"
                        if C[i] - decr > thr and C[i] <= thr:
                            Cnew_copy = Cnew.copy()
                            Cnew_copy[i] = C[i] - decr
                            Cnew_copy = 1 - Cnew_copy
                            update = simple_point_test(faces, i, Cnew_copy, thr)
                        # Or update the HMMF value if it is far from the threshold
                        else:
                            update = 1
                        # Update the HMMF and probability values
                        if update:
                            Cnew[i] = min([C[i] - decr, 1])
                            probs[i] = prob(wt_likelihood, L[i], wt_neighbors,
                                            Cnew[i], C[neighbors[i]])

            # Tally the number of vertices assigned probability values
            # and the number of HMMF values with probability > thr.
            # After iteration 1, compare current and previous values.
            # If the values are similar, increment end_flag.
            n_probs = sum(probs[indices])
            n_points = sum(C > thr)
            if count > 0:
                if n_points == n_points_previous:
                    diff_probs = (n_probs - n_probs_previous) / n_indices
                    if diff_probs < diff_threshold:
                        end_flag += 1

            # Reset for next iteration
            n_probs_previous = n_probs
            n_points_previous = n_points
            C = Cnew

            count += 1
            print(str(n_points) + ' vertices...')
    
        return C
