#!/usr/bin/python
"""
Compute distance measures.

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

#=================================
# Compute fundus likelihood values
#=================================
def compute_fundus_likelihood(mean_curvatures, depths, sulci, ind):
    """
    Compute fundus likelihood values.
    """
    threshold1 = .6
    threshold2 = .05
    threshold3 = .3
    high_map_value = .9
    depth_search_increment1 = 0.01
    depth_search_increment2 = 0.0001

    sulci(sulci != ind) = 0

    # Switch curvature values to opposite
    mean_curvatures = (-1 * mean_curvatures)

    depth_search = 0
    mass_left = 1

    # Find depth value where less than threshold1 are larger
    while(mass_left > threshold1):
        depth_search = depth_search + depth_search_increment1
        mass_left = sum(depths(sulci > 0) > depth_search) / sum(sulci > 0)
    half_depth = depth_search

    # Find depth value where less than threshold2 are larger
    while(mass_left > threshold2):
        depth_search = depth_search + depth_search_increment1
        mass_left = sum(depths(sulci > 0) > depth_search) / sum(sulci > 0)
    slope_depth = (-1/(depth_search - half_depth)) * log((1/high_map_value)-1)

    # Find slope similarly for curvature values
    mass_left = 1
    depth_search = 0
    while(mass_left > threshold3):
        depth_search = depth_search + depth_search_increment2
        mass_left = sum(mean_curvatures(sulci > 0) > depth_search) / sum(sulci > 0)
    slope_mean_curvatures = -1/depth_search * log((1/high_map_value) - 1)

    # Map curvature and depth values with sigmoidal functions to range [0,1]
    st2 = 1 ./ (1 + exp(-slope_mean_curvatures * mean_curvatures))
    st3 = 1 ./ (1 + exp(-slope_depth * (depths - half_depth)))

    L = st2 .* st3
    L(sulci ==0) = 0

    return L
