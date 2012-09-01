#!/usr/bin/python
"""
Find neighbors to a surface mesh vertex.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#---------------
# Find neighbors
#---------------
def find_neighbors(faces, index):
    """
    Find neighbors to a surface mesh vertex.

    For a set of surface mesh faces and the index of a surface vertex,
    find unique indices for neighboring vertices.

    Parameters
    ----------
    faces : surface mesh vertex indices [#faces x 3] numpy array
    index : index of surface vertex

    Returns
    -------
    N : list of indices of neighboring vertices

    """

    import numpy as np

    # Make sure argument is a numpy array
    if type(faces) != np.ndarray:
        faces = np.array(faces)

    # Create list of vertex indices sharing the same faces as "index"
    I = [faces[np.where(faces[:,i] == index)[0], :] for i in (0,1,2)]

    # Create single list from nested lists
    I = [int(x) for lst in I for sublst in lst for x in sublst]

    # Find unique indices not equal to "index"
    N = []; [N.append(x) for x in I if x not in N if x != index]

    return N


#===================
# Find anchor points
#===================
def find_anchors(vertices, L, min_directions, min_distance, thr):
    """
    Find anchor points.

    Assign maximum likelihood vertices as "anchor points"
    while ensuring that the anchor points are not close to one another.
    Anchor points are used to construct curves.

    Note: only the sulcus end points are strictly necessary
          (so that the fundus doesn't shrink)

    Parameters
    ----------
    vertices : [#vertices x 3] numpy array
    L : fundus likelihood values: [#vertices x 1] numpy array
    min_directions : [#vertices x 1] numpy array
    min_distance : minimum distance
    thr : likelihood threshold

    Returns
    -------
    anchors : list of subset of surface mesh vertex indices

    """

    import numpy as np
    from operator import itemgetter

    # Make sure arguments are numpy arrays
    if type(vertices) != np.ndarray:
        vertices = np.array(vertices)
    if type(L) != np.ndarray:
        L = np.array(L)
    if type(min_directions) != np.ndarray:
        min_directions = np.array(min_directions)

    max_distance = 2 * min_distance

    # Sort likelihood values and find indices for values above the threshold
    L_table = [[i,x] for i,x in enumerate(L)]
    L_table_sort = np.transpose(sorted(L_table, key=itemgetter(1)))[:, ::-1]
    IL = [int(L_table_sort[0,i]) for i,x in enumerate(L_table_sort[1,:])
          if x > thr]

    # Initialize anchors list with the index of the maximum likelihood value,
    # remove this value, and loop through the remaining high likelihoods
    anchors = [IL.pop(0)]
    for imax in IL:

        # Determine if there are any anchor points
        # near to the current maximum likelihood vertex
        i = 0
        found = 0
        while i < len(anchors) and found == 0:

            # Compute Euclidean distance between points
            D = np.linalg.norm(vertices[anchors[i], :] - vertices[imax, :])

            # If distance less than threshold, consider the point found
            if D < min_distance:
                found = 1
            # Compute directional distance between points if they are close
            elif D < max_distance:
                dirV = np.dot(vertices[anchors[i], :] - vertices[imax, :],
                              min_directions[anchors[i], :])
                # If distance less than threshold, consider the point found
                if np.linalg.norm(dirV) < min_distance:
                    found = 1

            i += 1

        # If there are no nearby anchor points,
        # assign the maximum likelihood vertex as an anchor point
        if not found:
            anchors.append(imax)

    return anchors
