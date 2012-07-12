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
    index: index of surface vertex

    Output:
    ------
    I: [#neighbors x 1] numpy array

    """
    # Create list of vertex indices sharing the same faces as "index"
    I = np.concatenate([x for x in faces if index in x])

    # Find unique indices not equal to "index"
    N = []
    [N.append(x) for x in I if x not in N if x != index]

    return N
