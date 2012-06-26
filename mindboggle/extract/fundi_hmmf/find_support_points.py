#!/usr/bin/python
"""
Compute distance measures.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

#-----------------------------
# Compute directional distance
#-----------------------------
def D = compDirectionalDist(P1,P2,U):
    """
    Compute directional distance.
    """
    dirV = np.dot(P1-P2,U)
    D = np.sqrt(sum(dirV.^2))

#---------------------------
# Compute Euclidean distance
#---------------------------
def d = compEucDist(A,B):
    """
    Compute Euclidean distance.
    """
    d = 0
    d = A-B
    d = d.^2
    d = np.sqrt(sum(d))

#====================
# Find support points
#====================
def P = findSupportPoints(vertices, values, Umin):
    """
    Find support points.
    """
    thr = 5

    P = zeros(size(values))
    checkList = P

    [maxVal maxInd] = max(values)
    checkList(maxInd) = 1
    P(maxInd) = 1
    values(maxInd) = -1

    while (min(checkList) == 0 && maxVal > .5):

        [maxVal maxInd] = max(values)
        checkList(maxInd) = 1
        values(maxInd) = -1

        currPVert = vertices(P>0,:)
        currU = Umin(:,P>0)
        i = 0
        found = 0
        while(i < size(currPVert,1) && found ==0):
            i = i + 1
            D = compEucDist(currPVert(i,:), vertices(maxInd,:))

            if (D >= thr && D<8):
                D = compDirectionalDist(currPVert(i,:), vertices(maxInd,:), currU(:,i))

            if (D < thr):
                found = 1

        if (found == 0):
            P(maxInd) = 1
