#!/usr/bin/python
"""
Connect the dots.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from simple_point_test import simple_point_test

#--------------------
# Compute probability
#--------------------
def prob(wl, li, q, qneighs, wt_neighbors):
    """
    Compute probability
    """
    p = -1 * (q * np.sqrt((wl-li)**2) + (wt_neighbors * sum((q - qneighs)**2)))
    return p

#-----------------------
# Test for simple points
#-----------------------
def simple_point_test(faces, ind, values):
    """
    Find anchor points

    Assign maximum likelihood vertices as "anchor points"
    for use in constructing fundus curves.
    Ensure tht anchor points are not close to one another.

    Inputs:
    ------
    vertices: [#vertices x 3] numpy array
    L: fundus likelihood values [#vertices x 1] numpy array
    min_distances: [#vertices x 1] numpy array
    distance_threshold: distance threshold
    max_distance: maximum distance

    Output:
    ------
    P: anchor points [#vertices x 1] numpy array

    """
    import numpy as np
    from find_neighbors import find_neigbhors

    thr = 0.5
    neighs = find_neighbors(faces,ind)
    neighvals = values(neighs)
    nOutside = sum(neighvals <= thr)
    nInside = sum(neighvals > thr)

    if (nOutside == 0 || nInside == 0)
        sp = 0
        return

    if (nOutside == 1)
        sp = 1
        return

    if (nInside == 1)
        sp = 1
        return

    neighsInside = neighs[neighvals > thr]
    nSets = len(neighsInside)

    sets = np.zeros(nSets,20)
    setLabels = np.transpose(range(1,nSets+1))

    values(ind) = -1

    for i in range(nSets):
        currNeighs = find_neighbors(faces,neighsInside[i])
        currNeighs = currNeighs(values(currNeighs)>thr)
        sets(i,1:size(currNeighs,1)) = np.transpose(currNeighs)
        sets(i,size(currNeighs,1) + 1) = neighsInside[i]

    change = 1
    while (change > 0):
        change = 0
        for i in range(nSets):
            for j in range(nSets):
               if (j>i && setLabels[i] != setLabels[j]):
                    iInds = np.unique(sets(i,:))
                    iInds = iInds(iInds > 0)
                    iSize = size(iInds,2)

                    jInds = np.unique(sets(j,:))
                    jInds = jInds(jInds > 0)
                    jSize = size(jInds,2)

                    totInds = unique([iInds jInds])
                    totSize = size(totInds,2)

                    if (totSize < iSize + jSize):

                        minLabel = min(setLabels[j],setLabels[i])

                        setLabels[i] = minLabel
                        setLabels[j] = minLabel

                        change = 1

    numberOfSeparateSets = size(unique(setLabels),1)
    if (numberOfSeparateSets < 2):
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
    C_threshold = 0.01  # ???? (positive)

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

            # ????
            mult = mult_init + count * mult_incr

            # Loop through vertices with positive likelihood values
            # and assign each a connectivity value (C[i])
            for i in I_positive:
                # Continue if connectivity value is greater than C_threshold
                if C[i] > C_threshold:
                    if any(indices == i):

                        # Find neighbors to this vertex
                        neighs = neighbors[indices == i, :]
                        neighs = np.unique(neighs[neighs > 0])

                        # Compute the probability for this point ????
                        q = max([C[i] - step_size, 0])
                        prob_down = prob(wl, L[i], q, C[neighs], wt_neighbors)
                        
                        # Compute a decrement...????
                        decr = mult * (prob_down - probs[i]) / step_size

                        # Test to update the connectivity value
                        # (default is not to update)
                        update = 0
                        # Update the connectivity value
                        # if it is far from L_threshold
                        if np.abs(C[i] - L_threshold) > np.abs(decr):
                            update = 1
                        # Otherwise...
                        elif decr > 0:
                            # If the dot value exceeds L_threshold
                            # do not update the connectivity value
                            if dots[i] > L_threshold:
                                update = 0
                            # Or if the dot value is within L_threshold
                            # and if the vertex is a "simple point"
                            # then update the connectivity value
                            else:
                                Cnew_copy = Cnew.copy()
                                Cnew_copy[i] = C[i] - decr
                                update = simple_point_test(faces, i, Cnew_copy)
                        elif decr < 0:
                            Cnew_copy = Cnew.copy()
                            Cnew_copy[i] = C[i] - decr
                            Cnew_copy = 1 - Cnew_copy
                            update = simple_point_test(faces, i, Cnew_copy)

                        # Update the connectivity and probability values
                        if update:
                            if decr > 0:
                                Cnew[i] = max([C[i] - decr, 0])
                            else:
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
                if diff_probs < diff_threshold and \
                   n_points == n_points_previous:
                    end_flag += 1

            # Reset for next iteration
            n_probs_previous = n_probs
            n_points_previous = n_points
            C = Cnew
            print('Iteration:', str(count) + ',', str(n_points), 'points')
    
        return C
