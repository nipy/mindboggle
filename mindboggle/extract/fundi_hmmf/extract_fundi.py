#!/usr/bin/python
"""
Extract fundi.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from find_neighbors import find_neigbhors
from extract_sulci import extract_sulci
from compute_fundus_likelihood import compute_fundus_likelihood

#=================
# Extract a fundus
#=================
def fundus = extract_fundus(L, Sulci, ind, vertices, faces, Umin):
    """
    Extract a single fundus.

    Determine minimum size for a sulcus from which we want to find fundi.
    This can be set to 0, the result will just have very short fundi.
    """

    minSulcusSize = 30

    L2 = np.zeros(size(L))
    L2(Sulci==ind) = L(Sulci==ind)

    sizeCheck = sum(L2>0.5)

    if (sizeCheck > minSulcusSize):
        Q = findSupportPoints(vertices, L2, Umin)

        # Compute list of neighbors
        inds = find(L2 > 0)
        m = len(inds,1)
        nList = np.zeros(m,15)
        for i in range(m):
            neighs = find_neighbors(faces, inds(i))
            nList(i, 1:len(neighs,1)) = np.transpose(neighs)
        shortInds = find(L2 > 0)

        # initialize all likelihood values within sulcus between 0.5 and 1.0. 
        # This is necessary to guarantee correct topology.
        initValues = (L2 + 1.001) / 2
        initValues(L2 == 0) = 0
        initValues(initValues > 1) = 1

        print(find(Q > 0.5))
        # Q = connectTheDotsV2(L2, initValues, vertices, faces, Q, nList, shortInds)
        if (sum(Q > 0.5) > 1):
            Q = connectTheDotsV3(L2, initValues, vertices, faces, Q, nList, shortInds)
            fundus = Q
        else:
            fundus = np.zeros(size(Q))
    else:
        print('Not Enough Values')
        print(sizeCheck)
        fundus = np.zeros(size(L2))

    return fundus

#==================
# Extract all fundi
#==================
def fundi = extract_all_fundi(vertices, faces, depths, mean_curvatures, min_directions):
    """
    Extract all fundi.

    Inputs:
    1. vertices has m x 3 elements
    2. faces has n x 3 elements
    3. depths: depth values in vector of size m x 1
    4. mean_curvatures: mean curvature values in vector of size m x 1
    5. min_directions: directions of minimum curvature in matrix of size 3 x m

    Parameters:
    -  depth_threshold: depth threshold for defining sulci
     
    Output:
     
    1. fundi: Matrix of size m x n_sulci, where n_sulci is the number of sulci.
              Each column represents a fundus for a specific sulcus. 
              Values range from 0 to 1; values above 0.5 are considered part of the fundus.
    """

    # Depth threshold for defining sulci
    depth_threshold = 0.2

    # Extract and topologically correct sulci
    sulci = extract_sulci(faces, depths, depth_threshold)
    n_sulci = np.max(sulci(:))  # number of sulci

    # Compute fundus likelihood values
    print('computing likelihood values')
    L = np.zeros(len(depths),1)
    for i in range(n_sulci):
        L = L + compute_fundus_likelihood(mean_curvatures, depths, sulci, i)

    # Extract individual fundi
    fundi = np.zeros(len(L), n_sulci)
    for i in range(n_sulci):
        fundi(:,i) = extract_fundus(L, sulci, i, vertices, faces, min_directions)

    return fundi
