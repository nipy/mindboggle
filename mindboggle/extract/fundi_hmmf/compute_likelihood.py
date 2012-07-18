#!/usr/bin/python
"""
Compute fundus likelihood values.


Authors:
    Yrjo Hame  .  yrjo.hame@gmail.com
    Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

verbose = 1  # 2 for debugging

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
    # Find depth and curvature values greater than
    # the values of some portion of the vertices
    depth_percentile1: percentile of depth values for computing "depth1"
    depth_percentile2: percentile of depth values for computing "depth2"
    curvature_percentile: percentile of curvature values for computing "curvature"
    high_map_value: used for computing the depth and curvature slope factors
    precision_limit: a large number used to prevent precision errors

    Output:
    ------
    likelihoods: likelihood values [#sulcus vertices x 1] numpy array

    Calls:
    -----
    percentile()

    """
    # Parameters
    depth_percentile1 = 0.4
    depth_percentile2 = 0.95
    curvature_percentile = 0.7
    high_map_value = 0.9
    precision_limit = 1000

    # Take the opposite of the curvature values
    curvatures = -curvatures

    # Find depth and curvature values greater than
    # the values of some portion of the vertices
    sort_depths = np.sort(depths)
    sort_curvatures = np.sort(curvatures)
    depth1 = percentile(sort_depths, depth_percentile1, key=lambda x:x)
    depth2 = percentile(sort_depths, depth_percentile2, key=lambda x:x)
    curvature = percentile(sort_curvatures, curvature_percentile, key=lambda x:x)
    if verbose == 2:
        print('    {0:.2f}, {1:.2f} depths greater than {2:.2f}, {3:.2f} of vertices'.
              format(depth1, depth2, depth_percentile1, depth_percentile2))
        print('    {0:.2f} curvature greater than {1:.2f} of vertices'.
              format(curvature, curvature_percentile))

    # Find slope for depth and curvature values
    # (while preventing precision errors for large slope values)
    slope_factor = np.log((1. / high_map_value) - 1)
    if depth2 == depth1:
        slope_depth = -precision_limit
        if verbose == 1:
            print('    (note: infinite slope depth)')
    else:
        slope_depth = -slope_factor / (depth2 - depth1)
        if slope_depth > precision_limit:
            slope_depth = precision_limit
            if verbose == 1:
                print('    (note: high +slope depth: ' + str(slope_depth) + ')')
        elif slope_depth < -precision_limit:
            slope_depth = -precision_limit
            if verbose == 1:
                print('    (note: high -slope depth: ' + str(slope_depth) + ')')

    if curvature == 0:
        slope_curvature = -precision_limit
        if verbose == 1:
            print('    (note: infinite slope curvature)')
    else:
        slope_curvature = -slope_factor / curvature
        if slope_curvature > precision_limit:
            slope_curvature = precision_limit
            if verbose == 1:
                print('    (note: high +slope curvature: ' + str(slope_curvature) + ')')
        elif slope_curvature < -precision_limit:
            slope_curvature = -precision_limit
            if verbose == 1:
                print('    (note: high -slope curvature: ' + str(slope_curvature) + ')')

    # Map values with sigmoidal function to range [0,1]
    st_depths = 1 / (1 + np.exp(-slope_depth * (depths - depth2)))
    st_curvatures = 1 / (1 + np.exp(-slope_curvature * curvatures))

    # Assign likelihood values to vertices
    likelihoods = st_depths * st_curvatures  # element-wise multiplication

    return likelihoods
