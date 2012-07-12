#!/usr/bin/python
"""
Compute fundus likelihood values.


Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com
    Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np
from scipy.stats import scoreatpercentile

print_debug = 1

#=================================
# Compute fundus likelihood values
#=================================
def compute_likelihood(depths, curvatures):
    """
    Compute fundus likelihood values.

    ????[Include an explanation of adaptive thresholding]

    Inputs:
    ------
    curvatures: mean curvature values [#sulcus vertices x 1] numpy array
    depths: depth values [#sulcus vertices x 1] numpy array

    Parameters:
    ----------
    Adaptive thresholding:
      depth_threshold1: ????
      depth_threshold2: ????
      curvature_threshold: ????
      high_map_value: ????
    Increments to reduce computation time:
      increment1
      increment2

    Output:
    ------
    likelihoods: likelihood values [#sulcus vertices x 1] numpy array

    """
    # Parameters for adaptive thresholding
    depth_threshold1 = 0.6
    depth_threshold2 = 0.05
    curvature_threshold = 0.3
    high_map_value = 0.9

    # Increments to reduce computation time
    depth_increment = 0.01
    curvature_increment = 0.0001

    slope_factor = np.log((1. / high_map_value) - 1)

    # Take the opposite of the curvature values
    curvatures = -curvatures

    # Find depths and curvature values
    len_object = np.float(len(depths))

    """
    # Alternative to resetting the threshold values below?
    depth_found = scoreatpercentile(depths, depth_threshold2)
    slope_depth = -slope_factor / (search - depth_found)
    # Map depth values with sigmoidal function to range [0,1]
    st_depths = 1 / (1 + np.exp(-slope_depth * (depths - depth_found)))
    """

    # Find depth value where less than depth_threshold1 of vertices are deeper
    mass_left = 1
    search = 0
    while mass_left > depth_threshold1:
        search += depth_increment
        mass_left = sum(depths > search) / len_object
    depth_found = search
    if print_debug:
        print(str(depth_threshold1) +
              ' of vertices are deeper than ' + str(search))

    # Find depth value where less than depth_threshold2 of sulcus vertices are deeper
    while mass_left > depth_threshold2:
        search += depth_increment
        mass_left = sum(depths > search) / len_object
    if search == depth_found:
        st_depths = np.zeros(len_object)
    else:
        slope_depth = -slope_factor / (search - depth_found)
        # Map depth values with sigmoidal function to range [0,1]
        st_depths = 1 / (1 + np.exp(-slope_depth * (depths - depth_found)))
    if print_debug:
        print(str(depth_threshold2) +
              ' of vertices are deeper than ' + str(search))

    # Find slope for curvature values
    mass_left = 1
    search = 0
    while mass_left > curvature_threshold:
        search += curvature_increment
        mass_left = sum(curvatures > search) / len_object
    if print_debug:
        print(str(curvature_threshold) +
              ' of vertices have greater curvature than ' + str(search))
    slope_curvature = -slope_factor / search
    # Prevent precsion errors
    if slope_curvature > 1000:
        if print_debug:
            print('(high slope curvature: ' + str(slope_curvature) + ')')
        curvatures[curvatures < 0] = np.Inf
        curvatures[curvatures > 0] = 0
    # Map curvature values with sigmoidal function to range [0,1]
    st_curvatures = 1 / (1 + np.exp(-slope_curvature * curvatures))

    # Assign likelihood values to vertices
    likelihoods = st_depths* st_curvatures  # element-wise multiplication

    return likelihoods
