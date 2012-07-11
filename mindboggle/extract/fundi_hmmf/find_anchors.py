#!/usr/bin/python
"""
Assign maximum likelihood vertices as "anchor points".

Anchor points are used to construct fundus curves.

Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com
    Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

#===================
# Find anchor points
#===================
def find_anchors(vertices, L, min_directions, thr, min_distance, max_distance):
    """
    Find anchor points.

    Assign maximum likelihood vertices as "anchor points"
    for use in constructing fundus curves.
    Ensure that the anchor points are not close to one another.

    ????
    Should perhaps change this code because of instability.
    Only the sulcus end points are necessary,
    so that the fundus doesn't shrink.

    Inputs:
    ------
    vertices: [#vertices x 3] numpy array
    L: fundus likelihood values: [#vertices x 1] numpy array
    min_directions: [#vertices x 1] numpy array
    thr: likelihood threshold
    min_distance: minimum distance
    max_distance: maximum distance

    Output:
    ------
    anchors: list of anchors (subset of surface mesh vertex indices)

    """

    # Find maximum likelihood value
    maxL = max(L)
    imax = np.where(L == maxL)[0][0]  #[i for i,x in enumerate(L) if x == maxL]

    # Initialize anchor point list and reset maximum likelihood value
    anchors = [imax]
    L[imax] = -1

    # Loop through high-likelihood vertices
    while maxL > thr:

        # Find and reset maximum likelihood value
        maxL = max(L)
        imax = np.where(L == maxL)[0][0]  #[i for i,x in enumerate(L) if x == maxL]
        L[imax] = -1

        # Find anchor points close to vertex with maximum likelihood value
        i = 0
        found = 0
        while i < len(anchors) and found == 0:

            # Compute Euclidean distance between points
            D = np.linalg.norm(vertices[anchors[i], :] - vertices[imax, :])

            # If distance less than threshold, consider the point found
            if D < min_distance:
                found = 1
            # Compute directional distance between points if they are close
            elif max_distance > D >= min_distance:
                dirV = np.dot(vertices[anchors[i], :] - vertices[imax, :],
                              min_directions[anchors[i], :])
                D = np.linalg.norm(dirV)
                # If distance less than threshold, consider the point found
                if D < min_distance:
                    found = 1

            i += 1

        # If there are no nearby anchor points,
        # assign the maximum likelihood vertex as an anchor point
        if not found:
            anchors.append(imax)

    return anchors
