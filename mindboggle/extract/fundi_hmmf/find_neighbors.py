#!/usr/bin/python
"""
Find neighbors to a surface mesh vertex.

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

