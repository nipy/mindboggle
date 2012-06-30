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
    I = [faces[np.where(faces[:,i] == index)[0]].tolist() for i in range(3)]

    # Create single list from nested lists
    I = [int(item) for sublist in I for subsublist in sublist for item in subsublist]

    if len(I) > 0:

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

    Calls:
    -----
    find_neighbors()

    """

    # Parameters
    min_sulcus_size = 50

    # Remove small sulci and renumber remaining sulcus vertices
    print('Number of sulci before pruning:', str(max(sulci)))
    counter = 0
    while counter < np.max(sulci):
        counter += 1
        if sum(sulci == counter) < min_sulcus_size:
            sulci[sulci == counter] = 0
            sulci[sulci > counter] = sulci[sulci > counter] - 1
            counter -= 1
    print('Number of sulci after pruning:', str(max(sulci)))

    # Initialize holes and hole seeds
    holes = np.zeros(len(sulci))
    seeds = np.where(sulci == 0)[0]
    n_seeds = len(seeds)

    # Loop until there are no more hole seeds
    counter = 0
    while n_seeds > 0:
        counter += 1

        # Create a new array the size of the sulci array
        # and mark a random (inner seed) element (with "2")
        TEMP = np.zeros(len(sulci))
        rseed = np.round(np.random.rand * (n_seeds - 1)) + 1
        TEMP[seeds[rseed]] = 2
        new_size = sum(TEMP > 1)

        # Grow seed region (hole) until no more connected points available
        while new_size > 0:

            # Identify and reset inner seeds
            inner_seeds = np.where(TEMP == 2)[0]
            TEMP[TEMP == 2] = 1

            # Mark neighbors to inner seeds
            for inner_seed in inner_seeds:
                neighs = find_neighbors(faces, inner_seed)
                neighs = neighs[sulci[neighs] == 0]
                neighs = neighs[TEMP[neighs] == 0]
                TEMP[neighs] = 2
            new_size = sum(TEMP > 1)

        # Assign counter number to hole
        holes[TEMP > 0] = counter

        # Assign 0.5 to hole vertices in "sulci" array to ignore in loop
        # and find new hole seeds
        sulci[holes > 0] = 0.5
        seeds = np.where(sulci == 0)[0]
        n_seeds = len(seeds)

        # Display current hole size
        print('Hole size:', str(sum(holes == counter)))

    # Remove hole values from sulci array (new ones added below)
    sulci[sulci < 1] = 0

    # Find the maximum hole size (the background) to ignore below
    max_hole_size = 0
    for i in range(max(holes)):
        if sum(holes == i) > max_hole_size:
            max_hole_size = sum(holes == i)

    # Look for vertices that have a sulcus label and are
    # connected to any of the vertices in the current hole,
    # and assign the hole the maximum label number
    for i in range(max(holes)):
        found = 0
        hole_indices = np.where(holes == i)[0]
        # Ignore the hole with the largest number of vertices (background)
        if len(hole_indices) < max_hole_size:
            # Loop until a labeled neighbor is found
            j = 0
            while not found:
                j += 1
                # Find neighboring vertices to the hole
                neighs = find_neighbors(faces, hole_indices[j])
                # If there are any neighbor sulcus labels,
                # assign the sulci hole vertices the maximum neighbor label
                # and end the while loop
                nIdx = max(sulci[neighs])
                if nIdx > 0:
                    sulci[hole_indices] = nIdx
                    found = 1
                    print('Found hole', str(i))

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
        rseed = int(np.round(np.random.rand() * (n_seeds - 1)))
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
                    neighbors = neighbors[depths[neighbors] > depth_threshold]
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
        print('Number of sulci: ' + str(counter))

    # Fill sulcus holes to preserve topology
    sulci = fill_sulcus_holes(faces, sulci)

    # Compute the number of sulci
    n_sulci = np.max(sulci)

    return sulci, n_sulci
