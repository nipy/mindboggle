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

#========================
# Segment surface patches
#========================
def segment_surface(faces, seeds, N, min_patch_size):
    """
    Segment a surface into contiguous patches (seed region growing).

    Inputs:
    ------
    faces: surface mesh vertex indices [#faces x 3]
    seeds: mesh vertex indices for vertices to be segmented [#seeds x 1]
    N: #vertices total
    min_patch_size: minimum size of segmented set of vertices

    Output:
    ------
    segments: label indices for patches: [#seeds x 1] numpy array
    max_patch: index for largest segmented set of vertices

    Calls:
    -----
    find_neighbors(): numpy array of indices

    """

    # Initialize sulci and seeds (indices of deep vertices)
    segments = np.zeros(N)
    n_seeds = len(seeds)

    # Remove faces that do not contain seeds to speed up computation
    fs = frozenset(seeds)
    faces_seeds = [lst for lst in faces if fs.intersection(lst)]
    faces_seeds = np.reshape(np.ravel(faces_seeds), (-1, 3))
    print('Reduced ' + str(len(faces)) + ' to ' +\
          str(len(faces_seeds)) + ' faces')

    # Loop until all seed vertices segmented
    print('Grow ' + str(n_seeds) + ' seed vertices...')
    max_patch_size = 0
    max_patch = 0
    counter = 0
    TEMP0 = np.zeros(N)
    while n_seeds > min_patch_size:
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
                neighbors = find_neighbors(faces_seeds, index)
                # Select neighbors that have not been previously selected
                if len(neighbors) > 0:
                    neighbors = neighbors[TEMP[neighbors] == 0]
                    TEMP[neighbors] = 2
                    if len(Inew) > 0:
                        Inew = np.concatenate((Inew, neighbors))
                    else:
                        Inew = neighbors
                    # Continue looping
                    loop = 1
            I = Inew

        # Find region grown from seed
        Ipatch = np.where(TEMP > 0)[0]

        # Disregard vertices already visited
        seeds = list(frozenset(seeds).difference(Ipatch))
        n_seeds = len(seeds)

        # Assign counter number to segmented patch
        # if patch size is greater than min_patch_size
        size_patch = len(Ipatch)
        if size_patch > min_patch_size:
            counter += 1
            segments[Ipatch] = counter
            # Display current number and size of patch
            print('Segmented patch ' + str(counter) + ': ' + str(size_patch) +\
                  ' vertices.  ' + str(n_seeds) + ' remaining...')

            # Find the maximum hole size (the background) to ignore below
            if size_patch > max_patch_size:
                max_patch_size = size_patch
                max_patch = counter

    return segments, max_patch

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

    # Remove faces that do not contain seeds to speed up computation
    fs = frozenset(seeds)
    faces_holes = [lst for lst in faces if fs.intersection(lst)]
    faces_holes = np.reshape(np.ravel(faces_holes), (-1, 3))
    print('Reduced ' + str(len(faces)) + ' to ' +\
          str(len(faces_holes)) + ' faces')

    # Loop until there are no more hole seeds
    max_hole_size = 0
    max_hole_count = 0
    counter = 0
    TEMP0 = np.zeros(len(sulci))
    while n_seeds > 0:

        # Create a new array the size of the sulci array
        # and mark a random (inner seed) element (with "2")
        TEMP = np.copy(TEMP0)
        #I = [seeds[round(np.random.rand() * (n_seeds - 1))]]
        I = [seeds[0]]

        # Grow region (hole) about the seed vertex until
        # there are no more connected seed vertices available.
        # Continue loop if there are newly selected neighbors.
        loop = 1
        while loop:
            loop = 0
            TEMP[I] = 1
            Inew = np.array([])
            # Find neighbors for each selected seed vertex
            for index in I:
                neighbors = find_neighbors(faces_holes, index)
                # Select neighbors that have not been previously selected
                if len(neighbors) > 0:
                    neighbors = neighbors[TEMP[neighbors] == 0]
                    TEMP[neighbors] = 2
                    if len(Inew) > 0:
                        Inew = np.concatenate((Inew, neighbors))
                    else:
                        Inew = neighbors
                    loop = 1
                # Continue looping
            I = Inew

        # Find hole region grown from seed
        Ihole = np.where(TEMP > 0)[0]

        # Disregard vertices already visited
        seeds = list(frozenset(seeds).difference(Ihole))
        n_seeds = len(seeds)

        # Assign counter number to hole
        size_hole = len(Ihole)
        counter += 1
        holes[Ihole] = counter
        # Display current number and sizes of holes
        if size_hole > 1:
            print('Hole ' + str(counter) + ': ' + str(size_hole) +\
                  ' vertices.  ' + str(n_seeds) + ' remaining...')

        # Find the maximum hole size (the background) to ignore below
        if size_hole > max_hole_size:
            max_hole_size = size_hole
            max_hole_count = counter

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
    sulci: label indices for sulci: [#vertices x 1] numpy array
    n_sulci:  #sulci [int]

    Calls:
    -----
    find_neighbors(): numpy array of indices
    fill_sulcus_holes()

    """

    # Initialize sulci and seeds (indices of deep vertices)
    N = len(depths)
    sulci = np.zeros(N)
    seeds = np.where(depths > depth_threshold)[0]
    n_seeds = len(seeds)

    # Remove faces that do not contain seeds to speed up computation
    fs = frozenset(seeds)
    faces_seeds = [lst for lst in faces if fs.intersection(lst)]
    faces_seeds = np.reshape(np.ravel(faces_seeds), (-1, 3))
    print('Reduced ' + str(len(faces)) + ' to ' +\
          str(len(faces_seeds)) + ' faces')

    # Loop until all seed vertices included in sulci
    print('Grow ' + str(n_seeds) + ' sulcus seed vertices...')
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
                neighbors = find_neighbors(faces_seeds, index)
                # Select neighbors that have not been previously selected
                if len(neighbors) > 0:
                    neighbors = neighbors[TEMP[neighbors] == 0]
                    TEMP[neighbors] = 2
                    if len(Inew) > 0:
                        Inew = np.concatenate((Inew, neighbors))
                    else:
                        Inew = neighbors
                    loop = 1
            # Continue looping
            I = Inew

        # Find sulcus region grown from seed
        Isulcus = np.where(TEMP > 0)[0]

        # Disregard vertices already visited
        seeds = list(frozenset(seeds).difference(Isulcus))
        n_seeds = len(seeds)

        # Assign counter number to sulcus
        # if sulcus size is greater than min_sulcus_size
        size_sulcus = len(Isulcus)
        if size_sulcus > min_sulcus_size:
            counter += 1
            sulci[Isulcus] = counter
            # Display current number and sizes of sulci
            print('Sulcus ' + str(counter) + ': ' + str(size_sulcus) +\
                  ' vertices.  ' + str(n_seeds) + ' remaining...')

    # Fill sulcus holes to preserve topology
    n_sulci = np.max(sulci)
    if n_sulci > 0:
        print('Filling holes in sulci...')
        sulci = fill_sulcus_holes(faces, sulci)
    else:
        import sys
        sys.exit('There are no sulci.')

    # Return sulci and the number of sulci
    return sulci, n_sulci

