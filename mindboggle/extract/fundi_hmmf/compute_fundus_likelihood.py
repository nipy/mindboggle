#!/usr/bin/python
"""
Compute fundus likelihood values.

Input


Authors:
Yrjo Hame  .  yrjo.hame@gmail.com  (original Matlab code)
Arno Klein  .  arno@mindboggle.info  (translated to Python)

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

#=================================
# Compute fundus likelihood values
#=================================
def compute_fundus_likelihood(curvatures, depths, sulci, sulcus_index):
    """
    Compute fundus likelihood values.

    ????

    Inputs:
    ------
    curvatures: mean curvature values [#vertices x 1] numpy array
    depths: depth values [#vertices x 1] numpy array
    sulci: sulcus values [#vertices x 1] numpy array
    sulcus_index: sulcus number [int]

    Parameters:
    ----------
    threshold1: ????
    threshold2: ????
    threshold3: ????
    high_map_value: ????
    increment1: ????
    increment2: ????

    Output:
    ------
    L: fundus likelihood values [#vertices x 1] numpy array

    """
    # must these be fixed????
    # ???? explain adaptive thresholding
    threshold1 = 0.6
    threshold2 = 0.05
    threshold3 = 0.3
    high_map_value = 0.9

    increment1 = 0.01
    increment2 = 0.0001

    slope_factor = np.log((1/high_map_value)-1)

    # Take the opposite of the curvature values
    curvatures = (-1 * curvatures)

    # Retain only sulcus numbers for the sulcus corresponding to sulcus_index
    sulci[sulci != sulcus_index] = 0
    Igt0 = sulci > 0
    # Find sulcus depths and curvature values
    len_sulcus = sum(Igt0)
    sulcus_depths = depths[Igt0]
    sulcus_curvatures = curvatures[Igt0]

    # Find depth value where less than threshold1 are larger
    mass_left = 1
    depth_search = 0
    while(mass_left > threshold1):
        depth_search = depth_search + increment1
        mass_left = sum(sulcus_depths > depth_search) / len_sulcus
    depth_found = depth_search

    # Find depth value where less than threshold2 are larger
    while(mass_left > threshold2):
        depth_search = depth_search + increment1
        mass_left = sum(sulcus_depths > depth_search) / len_sulcus
    slope_depth = -slope_factor / (depth_search - depth_found)

    # Find slope for curvature values
    mass_left = 1
    depth_search = 0
    while(mass_left > threshold3):
        depth_search = depth_search + increment2
        mass_left = sum(sulcus_curvatures > depth_search) / len_sulcus
    slope_curvatures = -slope_factor / depth_search

    # Map curvature and depth values with sigmoidal functions to range [0,1]
    st2 = 1 / (1 + np.exp(-slope_curvatures * curvatures))
    st3 = 1 / (1 + np.exp(-slope_depth * (depths - depth_found)))

    L = st2 * st3  # element-wise multiplication
    L[sulci == 0] = 0

    return L
