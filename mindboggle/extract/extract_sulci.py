#!/usr/bin/python
"""
Extract sulci from a VTK surface mesh and fill holes in them.

Inputs:
    faces: surface mesh vertex indices [n x 3]
    depths: depth values [m x 1]
    depth_threshold: depth threshold for defining sulci
     
Output:
    sulci: [n_sulci x 1], where n_sulci is the number of sulci

Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
    Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from find_neighbors import find_neighbors

#------------------
# Fill sulcus holes
#------------------
def fill_sulcus_holes(faces, sulci):
    """
    Fill sulcus holes.
    """
    print('Number of sulci before pruning: ' + str(max(sulci)))

    counter = 0

    while (counter < max(sulci)):
        counter = counter + 1
        ns = sum(sulci == counter)
        if (ns < 50):
            sulci(sulci == counter) = 0
            sulci(sulci > counter) = sulci(sulci > counter) - 1
            counter = counter - 1

    print('Number of sulci after pruning: ' str(max(sulci)))

    Holes = np.zeros(size(sulci))

    seeds = find(sulci == 0)

    seedSize = size(seeds,1)
    counter = 0

    while (seedSize > 0):
        counter = counter + 1
        TEMP = np.zeros(size(sulci))

        rseed = round(rand*(seedSize-1)) + 1

        TEMP(seeds(rseed,1),1) = 2
        newSize = size(find(TEMP>1))

        # grow region until no more connected points available
        while(newSize > 0):
            indsList = find(TEMP == 2)

            TEMP(TEMP == 2) = 1

            for i = 1:size(indsList,1):
                neighs = find_neighbors(faces,indsList(i,1))
                neighs = neighs(sulci(neighs) == 0)
                neighs = neighs(TEMP(neighs) == 0)
                TEMP(neighs) = 2

            newSize = size(find(TEMP>1))

        Holes(TEMP > 0) = counter

        sulci(Holes > 0) = .5
        seeds = find(sulci == 0)

        seedSize = size(seeds,1)

        # display current sulcus size
        print(size(find(Holes == counter),1))

    sulci(sulci < 1) = 0

    for i in range(max(Holes)):
        found = 0
        currInds = find(Holes == i)

        if (size(currInds,1) < 10000):
            j = 0
            while (found == 0):
                j = j + 1
                neighs = find_neighbors(faces, currInds(j,1))
                nIdx = max(sulci(neighs))
                if (nIdx > 0):
                    sulci(currInds) = nIdx
                    found = 1
                    print('Found ' + str(i))

    return sulci

#==============
# Extract sulci
#==============
def extract_sulci(faces, depths, depth_threshold=0.2):
    """
    Extract sulci.

    Inputs:
    faces: surface mesh vertex indices [n x 3]
    depths: depth values [m x 1]
    depth_threshold: depth threshold for defining sulci
     
    Output:
    sulci: [n_sulci x 1], where n_sulci is the number of sulci

    """
    import numpy as np
    from find_neighbors import find_neigbhors

    sulci = np.zeros(size(depths))

    seeds = find(depths > depth_threshold)

    seed_size = size(seeds,1)
    counter = 0

    # loop until all points included in sulci
    while (seed_size > 0):
        counter = counter + 1
        TEMP = np.zeros(size(depths))

        #choose random seed point (selection does not affect result)
        rseed = np.round(rand*(seed_size-1)) + 1

        TEMP(seeds(rseed,1),1) = 2
        new_size = size(find(TEMP>1))

        # grow region until no more connected points available
        while(new_size > 0):
            indices = find(TEMP == 2)

            TEMP(TEMP == 2) = 1

            for i in range(size(indices,1)):
                neighbors = find_neighbors(faces, indices(i,1))
                neighbors = neighbors(depths(neighbors) > depth_threshold)
                neighbors = neighbors(TEMP(neighbors) == 0)
                TEMP(neighbors) = 2
            new_size = size(find(TEMP>1))
        sulci(TEMP > 0) = counter

        depths(sulci > 0) = depth_threshold - .5
        seeds = find(depths > depth_threshold)

        seed_size = size(seeds,1)

        # display current sulcus size
        print(size(find(sulci == counter),1))

    sulci = fill_sulcus_holes(faces, sulci)

    return sulci
