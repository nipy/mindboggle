#!/usr/bin/python
"""
Extract sulci from a VTK surface mesh and fill holes in them.

Inputs:
    faces: surface mesh vertex indices [#faces x 3]
    depths: depth values [#vertices x 1]
    depth_threshold: depth threshold for defining sulci

Parameters:
    depth_increment: depth increment for assigning vertices to sulci

Output:
    sulci: sulcus numbers [#vertices x 1]
    n_sulci:  #sulci


Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

debug_verbose = 1

#---------------
# Find neighbors
#---------------
def find_neighbors(faces, index):
    """
    For a set of surface mesh faces and the index of a surface vertex,
    find unique indices for neighboring vertices.

    Inputs:
    ------
    faces: surface mesh vertex indices [#faces x 3] numpy array
    index: index of surface vertex

    Output:
    ------
    I: [#neighbors x 1] numpy array

    """
    # Create list of vertex indices sharing the same faces as "index"
    I = [faces[np.where(faces[:,i] == index)[0]].tolist() for i in range(3)]

    # Create single list from nested lists
    I = [int(item) for sublist in I for subsublist in sublist for item in subsublist]

    if len(I) > 0:

        # Find unique indices not equal to "index"
        I = np.unique(I)
        I = I[I != index]

    return I

#------------------
# Fill sulcus holes
#------------------
def fill_sulcus_holes(faces, sulci):
    """
    Fill sulcus holes.

    ????[include explanation of sulci within holes problem]

    Inputs:
    ------
    faces: surface mesh vertex indices [#faces x 3] numpy array
    sulci: [#vertices x 1] numpy array

    Output:
    ------
    sulci: [#vertices x 1] numpy array

    Calls:
    -----
    find_neighbors()

    """

    # Initialize holes and hole seeds
    holes = np.zeros(len(sulci))
    seeds = np.where(sulci == 0)[0]
    n_seeds = len(seeds)

    # Loop until there are no more hole seeds
    max_hole_size = 0
    max_hole_count = 0
    counter = 0
    while n_seeds > 0:
        counter += 1

        if debug_verbose:
            print('#counter: '+str(counter))
            print('#seeds: '+str(n_seeds))

        # Create a new array the size of the sulci array
        # and mark a random (inner seed) element (with "2")
        TEMP = np.zeros(len(sulci))
        rseed = round(np.random.rand() * n_seeds) - 1
        TEMP[seeds[rseed]] = 2
        new_size = sum(TEMP > 1)

        # Grow seed region (hole) until no more connected points available
        while new_size > 0:

            # Identify and reset inner seeds
            inner_seeds = np.where(TEMP == 2)[0]
            TEMP[TEMP == 2] = 1

            if debug_verbose:
                print('#inner_seeds: '+str(len(inner_seeds)))

            # Mark neighbors to inner seeds
            for inner_seed in inner_seeds:
                neighs = find_neighbors(faces, inner_seed)
                if len(neighs) > 0:
                    neighs = neighs[sulci[neighs] == 0]
                    if len(neighs) > 0:
                        neighs = neighs[TEMP[neighs] == 0]
                        if len(neighs) > 0:
                            TEMP[neighs] = 2
            new_size = sum(TEMP > 1)

        # Assign counter number to hole
        hole_bool = TEMP > 0
        holes[hole_bool] = counter
        # Find the maximum hole size (the background) to ignore below
        hole_size = sum(hole_bool)
        if hole_size > max_hole_size:
            max_hole_size = hole_size
            max_hole_count = counter

        # Assign 0.5 to hole vertices in "sulci" array to ignore in loop
        # and find new hole seeds
        sulci[hole_bool] = 0.5
        seeds = np.where(sulci == 0)[0]
        n_seeds = len(seeds)

    # Remove hole values from sulci array (new ones added below)
    sulci[sulci < 1] = 0

    # Ignore the largest hole (the background)
    if max_hole_size > 0:
        holes[holes == max_hole_count] = 0

    # Look for vertices that have a sulcus label and are
    # connected to any of the vertices in the current hole,
    # and assign the hole the maximum label number
    for i in np.unique(holes):
        if i > 0:
            found = 0
            hole_indices = np.where(holes == i)[0]
            # Loop until a labeled neighbor is found
            for j in range(len(hole_indices)):
                # Find neighboring vertices to the hole
                neighs = find_neighbors(faces, hole_indices[j])
                # If there are any neighbor sulcus labels,
                # assign the sulci hole vertices the maximum neighbor label
                # and end the while loop
                if any(neighs):
                    nIdx = max(sulci[neighs])
                    if nIdx > 0:
                        sulci[hole_indices] = nIdx
                        break

    return sulci

#==============
# Extract sulci
#==============
def extract_sulci(faces, depths, depth_threshold=0.2, min_sulcus_size=50):
    """
    Extract sulci.

    Inputs:
    ------
    faces: surface mesh vertex indices [#faces x 3]
    depths: depth values [#vertices x 1]
    depth_threshold: depth threshold for defining sulci
    min_sulcus_size: minimum sulcus size

    Output:
    ------
    sulci: [#vertices x 1] numpy array
    n_sulci:  #sulci [int]

    Calls:
    -----
    find_neighbors(): numpy array of indices
    fill_sulcus_holes()

    """

    remove_faces = 1

    # Initialize sulci and seeds (indices of deep vertices)
    N = len(depths)
    sulci = np.zeros(N)
    seeds = np.where(depths > depth_threshold)[0]
    n_seeds = len(seeds)

    # Remove faces that do not contain seeds to speed up computation
    if remove_faces:
        faces_sulci = faces[~np.all((faces < min(seeds)) + (faces > max(seeds)), axis=1)]
#
#        faces_sulci = np.reshape(np.ravel([lst for lst in faces if len(np.intersect1d(seeds, lst)) > 0]), (-1,3))
#
#        print('Reduced ' + str(len(faces)) + ' to ' + \
#              str(len(faces_sulci)) + ' faces')


#        faces_sulci = faces.copy()
#        for seed in seeds:
#            faces_sulci *= (faces - seed)
#        faces_sulci = faces[np.sum(faces_sulci == 0, axis = 1) > 0]

        #boolsum = faces_sulci == seeds[0]
        #for seed in seeds:
        #    boolsum += faces_sulci == seed
        #faces_sulci = faces_sulci[np.sum(boolsum, axis = 1)]
        print('Reduced ' + str(len(faces)) + ' to ' + \
              str(len(faces_sulci)) + ' faces')

    # Loop until all seed vertices included in sulci
    print(str(n_seeds) + ' sulcus seed vertices to grow...')
    counter = 0
    TEMP0 = np.zeros(N)
    while n_seeds > min_sulcus_size:
        TEMP = np.copy(TEMP0)

        # Select a seed vertex (selection does not affect result)
        #I = [seeds[round(np.random.rand() * (n_seeds - 1))]]
        I = [seeds[0]]

        # Grow region about the seed vertex until 
        # there are no more connected seed vertices available.
        # Continue loop if there are newly selected neighbors.
        loop = 1
        while loop:
            loop = 0
            TEMP[I] = 1
            Inew = np.array([])
            # Find neighbors for each selected seed vertex
            for index in I:
                neighbors = find_neighbors(faces_sulci, index)
                # Select neighbors that have not been previously selected
                if len(neighbors) > 0:
                    neighbors = neighbors[TEMP[neighbors] == 0]
                    # Select neighbors deeper than the depth threshold
                    if len(neighbors) > 0:
                        neighbors = neighbors[depths[neighbors] > depth_threshold]
                        TEMP[neighbors] = 2
                        if len(Inew) == 0:
                            Inew = neighbors
                        else:
                            Inew = np.concatenate((Inew, neighbors))
                        loop = 1
            # Continue looping
            I = Inew

        # Find sulcus region grown from seed
        sulcus_bool = TEMP > 0

        # Disregard vertices already assigned values in sulci
        depths[sulcus_bool] = 0
        seeds = np.where(depths > depth_threshold)[0]
        n_seeds = len(seeds)

        # Assign the sulcus the counter number
        # if sulcus size is greater than min_sulcus_size
        size_sulcus = sum(sulcus_bool)
        if size_sulcus > min_sulcus_size:
            counter += 1
            sulci[sulcus_bool] = counter
            # Display current number and sizes of sulci
            print('Sulcus ' + str(counter) + ': ' + str(size_sulcus) +\
                  ' vertices.  ' + str(n_seeds) + ' remaining...')
        # Otherwise, assign small sulci indices 0.5 for later removal
        else:
            print('Candidate sulcus too small (' + str(size_sulcus) + ')')
            sulci[sulcus_bool] = 0.5

    # Remove small sulci (sulci indices assigned a value of 0.5)
    sulci[sulci == 0.5] = 0

    # Fill sulcus holes to preserve topology
    n_sulci = np.max(sulci)
    if n_sulci > 0:
        print('Filling holes in sulci...')
        sulci = fill_sulcus_holes(faces, sulci)
        #n_sulci = np.max(sulci)
    else:
        import sys
        sys.exit('There are no sulci.')

    # Return sulci and the number of sulci
    return sulci, n_sulci

