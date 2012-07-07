#!/usr/bin/python
"""
Compute fundus likelihood values.


Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com
    Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

#=================================
# Compute fundus likelihood values
#=================================
def compute_likelihood(sulci, sulcus_index, depths, curvatures):
    """
    Compute fundus likelihood values.

    ????[Include an explanation of adaptive thresholding]

    Inputs:
    ------
    curvatures: mean curvature values [#vertices x 1] numpy array
    depths: depth values [#vertices x 1] numpy array
    sulci: sulcus values [#vertices x 1] numpy array
    sulcus_index: sulcus number [int]

    Parameters:
    ----------
    Adaptive thresholding:
      threshold1: ????
      threshold2: ????
      threshold3: ????
      high_map_value: ????
    Increments to reduce computation time:
      increment1
      increment2

    Output:
    ------
    L: fundus likelihood values [#vertices x 1] numpy array

    """
    # Parameters for adaptive thresholding
    threshold1 = 0.6
    threshold2 = 0.05
    threshold3 = 0.3
    high_map_value = 0.9
    # Increments to reduce computation time
    increment1 = 0.01
    increment2 = 0.0001

    slope_factor = np.log((1. / high_map_value) - 1)

    # Take the opposite of the curvature values
    curvatures = -curvatures

    # Retain only sulcus numbers for the sulcus corresponding to sulcus_index
    sulcus_bool = sulci == sulcus_index

    # Find sulcus depths and curvature values
    len_sulcus = sum(sulcus_bool)
    sulcus_depths = depths[sulcus_bool]
    sulcus_curvatures = curvatures[sulcus_bool]

    # Find depth value where less than threshold1 are larger
    mass_left = 1
    depth_search = 0
    while mass_left > threshold1:
        depth_search += increment1
        mass_left = sum(sulcus_depths > depth_search) / len_sulcus
    depth_found = depth_search

    # Find depth value where less than threshold2 are larger
    while mass_left > threshold2:
        depth_search += increment1
        mass_left = sum(sulcus_depths > depth_search) / len_sulcus
    if depth_search == depth_found:
        slope_depth = -np.Inf
    else:
        slope_depth = -slope_factor / (depth_search - depth_found)

    # Find slope for curvature values
    mass_left = 1
    depth_search = 0
    while mass_left > threshold3:
        depth_search += increment2
        mass_left = sum(sulcus_curvatures > depth_search) / len_sulcus
    if not depth_search:
        slope_curvatures = -np.Inf
    else:
        slope_curvatures = -slope_factor / depth_search

    # Map curvature and depth values with sigmoidal functions to range [0,1]
    st2 = 1 / (1 + np.exp(-slope_curvatures * curvatures))
    st3 = 1 / (1 + np.exp(-slope_depth * (depths - depth_found)))

    L = st2 * st3  # element-wise multiplication
    L[~sulcus_bool] = 0

    return L
