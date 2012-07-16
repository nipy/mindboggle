#!/usr/bin/python
"""
Compute fundus likelihood values.


Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com
    Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

print_debug = 1

#============================
# Compute score at percentile
#============================
# http://code.activestate.com/recipes/511478/ (r2)
# Alternative scipy implementation:
# from scipy.stats import scoreatpercentile
# depth_found = scoreatpercentile(depths, depth_threshold2)

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if len(N) == 0:
        return None

    k = (len(N)-1) * percent
    f = np.floor(k)
    c = np.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)

    return d0 + d1


#=================================
# Compute fundus likelihood values
#=================================
def compute_likelihood(depths, curvatures):
    """
    Compute fundus likelihood values.

    Inputs:
    ------
    curvatures: mean curvature values [#sulcus vertices x 1] numpy array
    depths: depth values [#sulcus vertices x 1] numpy array

    Parameters:
    ----------
    Adaptive thresholding:
      depth_percentile1: ????
      depth_percentile2: ????
      curvature_percentile: ????
      high_map_value: ????

    Output:
    ------
    likelihoods: likelihood values [#sulcus vertices x 1] numpy array

    """
    # Parameters
    depth_percentile1 = 0.4
    depth_percentile2 = 0.95
    curvature_percentile = 0.7
    high_map_value = 0.9

    # Take the opposite of the curvature values
    curvatures = -curvatures

    # Find depth values where less than some threshold of vertices are deeper
    sort_depths = np.sort(depths)
    sort_curvatures = np.sort(curvatures)
    depth1 = percentile(sort_depths, depth_percentile1, key=lambda x:x)
    depth2 = percentile(sort_depths, depth_percentile2, key=lambda x:x)
    curvature = percentile(sort_curvatures, curvature_percentile, key=lambda x:x)
    if print_debug:
        print('    {0:.2f}, {1:.2f} of vertices deeper than {2:.2f}, {3:.2f} depths'.
              format(depth_percentile1, depth_percentile2, depth1, depth2))
        print('    {0:.2f} of vertices have greater curvature than {1:.2f}'.
              format(curvature_percentile, curvature))

    # Find slope for depth and curvature values
    slope_factor = np.log((1. / high_map_value) - 1)
    slope_depth = -slope_factor / (depth2 - depth1)
    slope_curvature = -slope_factor / curvature

    # Prevent precision errors
    if slope_curvature > 1000:
        if print_debug:
            print('    (high slope curvature: ' + str(slope_curvature) + ')')
        curvatures[curvatures < 0] = np.Inf
        curvatures[curvatures > 0] = 0

    # Map depth and curvature values with sigmoidal function to range [0,1]
    st_depths = 1 / (1 + np.exp(-slope_depth * (depths - depth2)))
    st_curvatures = 1 / (1 + np.exp(-slope_curvature * curvatures))

    # Assign likelihood values to vertices
    likelihoods = st_depths * st_curvatures  # element-wise multiplication

    return likelihoods
