#!/usr/bin/python
"""
Find neighbors

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

#===============
# Find neighbors
#===============
def inds = find_neighbors(faces,pointOfInterestIndex):
    """
    Find neighbors
    """

    inds = [faces(inds1,:);faces(inds2,:);faces(inds3,:)]
    inds = inds(:)
    inds = unique(inds)

    inds(inds == pointOfInterestIndex) = 0
    inds = inds(inds > 0)
