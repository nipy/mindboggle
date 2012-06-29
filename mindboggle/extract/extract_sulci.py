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
    sulci: [#vertices x 1] numpy array

    Output:
    ------
    I: [#vertices x 1] numpy array

    """
    # Create list of vertex indices sharing the same faces as "index"
    I = [faces[np.where(faces[:,i] == index)[0]][0].tolist() for i in range(3) \
         if len(np.where(faces[:,i] == index)[0]) > 0]

    # Create single list from nested lists
    I = [int(item) for sublist in I for item in sublist]

    # Find unique indices not equal to "index"
    I = np.unique(I)
    I[I != index]

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

    Parameters:
    ----------
    min_sulcus_size
    max_indices

    Calls:
    -----
    find_neighbors()

    """

    # Parameter
    min_sulcus_size = 50
    max_indices = 10000

    # Remove small sulci
    print('Number of sulci before pruning:', str(max(sulci)))
    counter = 0
    while counter < np.max(sulci):
        counter += 1
        ns = sum(sulci == counter)
        if ns < min_sulcus_size:
            sulci[sulci == counter] = 0
            sulci[sulci > counter] = sulci[sulci > counter] - 1
            counter -= 1
    print('Number of sulci after pruning: ' str(max(sulci)))

    # Initialize holes
    holes = np.zeros(len(sulci))
    seeds = np.where(sulci == 0)
    seed_size = len(seeds)

    counter = 0
    while seed_size > 0:
        counter = counter + 1
        TEMP = np.zeros(len(sulci))

        rseed = np.round(np.random.rand * (seed_size - 1)) + 1

        TEMP[seeds[rseed]] = 2
        new_size = sum(TEMP > 1)

        # grow region until no more connected points available
        while new_size > 0:
            indices = np.where(TEMP == 2)

            TEMP[TEMP == 2] = 1

            for index in indices:
                neighs = find_neighbors(faces, index)
                neighs = neighs[sulci[neighs] == 0]
                neighs = neighs[TEMP[neighs] == 0]
                TEMP[neighs] = 2

            new_size = sum(TEMP > 1)

        holes[TEMP > 0] = counter

        sulci[holes > 0] = 0.5
        seeds = np.where(sulci == 0)
        seed_size = len(seeds)

        # display current sulcus size
        print(sum(holes == counter))

    sulci[sulci < 1] = 0

    for i in range(max(holes)):
        found = 0
        current_indices = np.where(holes == i)

        if (len(current_indices) < max_indices):
            j = 0
            while found == 0:
                j += 1
                neighs = find_neighbors(faces, current_indices[j])
                nIdx = np.max(sulci[neighs])
                if nIdx > 0:
                    sulci[current_indices] = nIdx
                    found = 1
                    print('Found', str(i))

    return sulci

#==============
# Extract sulci
#==============
def extract_sulci(faces, depths, depth_threshold=0.2):
    """
    Extract sulci.

    Inputs:
    ------
    faces: surface mesh vertex indices [#faces x 3]
    depths: depth values [#vertices x 1]
    depth_threshold: depth threshold for defining sulci
     
    Parameters:
    ----------
    depth_increment: depth increment for assigning vertices to sulci
    
    Output:
    ------
    sulci: [#vertices x 1] numpy array
    n_sulci:  #sulci [int]

    Calls:
    -----
    find_neighbors(): numpy array of indices
    fill_sulcus_holes()

    """

    # Dummy input:
    #faces = np.round(10*np.random.rand(5,3)).astype(int)
    #depths = np.random.rand(10,1)

    depth_increment = 0.5 

    # Initialize sulci and seeds (indices of deep vertices)
    sulci = np.zeros(len(depths))
    seeds = np.where(depths > depth_threshold)[0]
    n_seeds = len(seeds)

    # Loop until all seed vertices included in sulci
    counter = 0
    while (n_seeds > 0):
        counter += 1

        # Select a random seed vertex (selection does not affect result)
        rseed = np.round(np.random.rand() * (n_seeds - 1))
        TEMP = np.zeros(len(depths))
        TEMP[seeds[rseed]] = 2

        # Grow region about the seed vertex until 
        # there are no more connected seed vertices available.
        # Continue loop if there are newly selected neighbors.
        while(sum(TEMP > 1) > 0):
            indices = np.where(TEMP == 2)[0]
            TEMP[indices] = 1
            # For each previously selected seed vertex
            for index in indices:
                # Find all neighbors deeper than the depth threshold
                neighbors = find_neighbors(faces, index)  # numpy array
                if any(neighbors):
                    neighbors = neighbors[depths[neighbors][0] > depth_threshold]
                    # Select the neighbors that have not been previously selected
                    if any(neighbors):
                        neighbors = neighbors[TEMP[neighbors] == 0]
                        TEMP[neighbors] = 2

        # Assign the seed region vertices the loop counter number
        sulci[TEMP > 0] = counter

        # Disregard vertices already assigned to sulci
        depths[sulci > 0] = depth_threshold - depth_increment
        seeds = np.where(depths > depth_threshold)[0]
        n_seeds = len(seeds)

        # Display current number of sulci
        print('Number of sulci:', str(counter))

    # Fill sulcus holes to preserve topology
    sulci = fill_sulcus_holes(faces, sulci)

    # Compute the number of sulci
    n_sulci = np.max(sulci)

    return sulci, n_sulci
