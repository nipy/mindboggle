#!/usr/bin/python
"""
Connect the dots.

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
def prob(wl, li, q, qneighs, wt_neighbors):
    """
    Compute probability????

    ????
    p = -q_i * sqrt(w-l_i)**2) - lambda * sum_j(q_i - q_j)**2


    term 1 promotes high likelihood values
    term 2 promotes smoothness of the HMMF values

    Inputs:
    ------
    wl: ????
    li: ????
    q: ????
    qneighs: ????
    wt_neighbors: ????

    Output:
    ------
    p: probability????

    """
    p = -1 * (q * np.sqrt((wl-li)**2) + (wt_neighbors * sum((q - qneighs)**2)))
    return p

#-----------------------
# Test for simple points
#-----------------------
def simple_point_test(faces, index, values, thr=0.5):
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

    Output:
    ------
    sp: simple point or not? (1 or 0)

    Calls:
    -----
    find_neighbors()

    """

    max_set_size = 20

    # Find neighbors to the input vertex and
    # count the number of "inside" and "outside" neighbors
    # (n_sets for neighbors greater than the likelihood threshold, thr,
    #  n_outside for neighbors less than or equal to thr)
    neighs = find_neighbors(faces, index)
    neighvals = values(neighs)
    n_sets = sum(neighvals > thr)
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
        neighs_inside = neighs[neighvals > thr]

        # Initialize "sets" numpy array [#neighbors>thr x max_set_size]
        sets = np.zeros((n_sets, max_set_size))
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
            if len_current > max_set_size - 1:
                len_current = max_set_size - 1
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
def connect_the_dots(L, L_init, faces, dots, neighbors, indices):
    """
    Connect dots (vertices in a surface mesh) to create curves.

    ???? explanation with HMMF description
    -- include equations and references parameters below as vars in the equation

    Inputs:
    ------
    L: likelihood values for a subset of vertices (zero for the rest):
       [#vertices x 1] numpy array
    L_init: initial likelihood values from 0 to 1: [#vertices x 1] numpy array
    vertices: [#vertices x 3] numpy array
    faces: [#faces x 3] numpy array
    dots: [#vertices x 1] numpy array
    indices: indices of subset of vertices: [#vertices_subset x 1] numpy array
    neighbors: indices of neighbors to each vertex in vertices_subset:
               [#vertices_subset x #max_neighbors] numpy array

    Parameters:
    ----------
    L_threshold
    step_size
    wt_neighbors
    wl
    mult_init

    Output:
    ------
    fundus: [#vertices x 1] numpy array

    Calls:
    -----
    prob()
    simple_point_test()

    """

    # Set parameters
    L_threshold = 0.5  # ????
    step_size = 0.05  # ????
    wt_neighbors = 0.4  # ????
    wl = 1.1  # ????
    mult_init = 0.02  # ????
    n_tries_no_change = 3  # ????
    max_count = 350  # ????
    diff_threshold = 0.0001  # ????
    mult_incr = 0.001  # ????
    C_threshold = 0.01  # ???? ...to fix at very low values, to make the optimization faster)

    # Initialize array of connected dots,
    # containing likelihood values greater than L_threshold,
    # and 1 for dots greater than L_threshold
    n_vertices = len(L)
    C = np.zeros(n_vertices)
    C[L_init > L_threshold] = L_init[L_init > L_threshold]
    C[dots > L_threshold] = 1

    # Continue if there are at least two candidate dots
    if sum(C > 0) >= 2:
        print('Initial candidates: ', str(sum(C > 0)))

        # Initialize new arrays of connected points and probability values
        Cnew = C.copy()
        probs = np.zeros(n_vertices)
        L_positive = L > 0
        n_positive = len(L_positive)
        I_positive = np.where(L_positive)
    
        # Assign probability values to each neighboring vertex
        for i in range(n_vertices):
            if any(indices == i):
                neighs = neighbors[indices == i, :]
                neighs = np.unique(neighs[neighs > 0])
                probs[i] = prob(wl, L[i], C[i], C[neighs], wt_neighbors)

        # Loop until count reaches max_count or until end_flag equals zero
        # (end_flag is used to allow the loop to continue even if there is
        #  no change for n_tries_no_change times)
        count = 0
        end_flag = 0
        while end_flag < n_tries_no_change and count < max_count:
            count += 1

            # Define a factor to multiply the probability gradient that will
            # increase the time-step size toward the end of the optimization
            mult = mult_init + count * mult_incr

            # Loop through vertices with positive likelihood values
            # and assign each a locally optimal HMMF value (C[i]) 
            # ????fuzzy continuous labeling????
            for i in I_positive:
                # Continue if HMMF value is greater than C_threshold
                # (to fix at very low values, to make the optimization faster)
                if C[i] > C_threshold:
                    if any(indices == i):

                        # Find neighbors to this vertex
                        neighs = neighbors[indices == i, :]
                        neighs = np.unique(neighs[neighs > 0])

                        # Compute the (negative) probability gradient for the HMMF value
                        # ????
                        q = max([C[i] - step_size, 0])
                        prob_decr = prob(wl, L[i], q, C[neighs], wt_neighbors)
                        decr = mult * (prob_decr - probs[i]) / step_size

                        # Test to update the HMMF value:
                        # Update the HMMF value
                        # if it is far from L_threshold
                        # ???? --- HAS TO CROSS THRESHOLD!!!
                        if np.abs(C[i] - L_threshold) > np.abs(decr):
                            update = 1
                        # Otherwise, if the decrement is positive,
                        # the dot value is less than or equal to L_threshold,
                        # and the vertex is a "simple point,"
                        # then update the HMMF value
                        elif decr > 0 and dots[i] <= L_threshold:
                            Cnew_copy = Cnew.copy()
                            Cnew_copy[i] = C[i] - decr
                            update = simple_point_test(faces, i, Cnew_copy, L_threshold)
                        # Or if the decrement is negative,
                        # and the vertex is a "simple point,"
                        # then update the HMMF value
                        elif decr < 0:
                            Cnew_copy = Cnew.copy()
                            Cnew_copy[i] = C[i] - decr
                            Cnew_copy = 1 - Cnew_copy
                            update = simple_point_test(faces, i, Cnew_copy, L_threshold)
                        # The default is not to update
                        else:
                            update = 0

                        # Update the connectivity and probability values
                        if update:
                            if decr > 0:
                                Cnew[i] = max([C[i] - decr, 0])
                            elif decr < 0:
                                Cnew[i] = min([C[i] - decr, 1])
                            probs[i] = prob(wl, L[i], Cnew[i], C[neighs], wt_neighbors)

            # Tally the number of vertices assigned probability values
            # and the number of points in C with probability > L_threshold.
            # After iteration 1, compare current and previous values.
            # If the values are similar, increment end_flag.
            n_probs = sum(probs[L_positive])
            if count > 0:
                n_points = sum(C > L_threshold)
                diff_probs = (n_probs - n_probs_previous) / n_positive
                if diff_probs < diff_threshold and n_points == n_points_previous:
                    end_flag += 1

            # Reset for next iteration
            n_probs_previous = n_probs
            n_points_previous = n_points
            C = Cnew
            print('Iteration:', str(count) + ',', str(n_points), 'points')
    
        return C
