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
    depth_percentile_low: percentile of depth values for computing "depth_percentile_low"
    depth_percentile_high: percentile of depth values for computing "depth_percentile_high"
    curvature_percentile: percentile of curvature values for computing "curvature"
    high_map_value: used for computing the depth and curvature slope factors
    precision_limit: a large number in place of infinity to prevent precision errors

    Output:
    ------
    likelihoods: likelihood values [#sulcus vertices x 1] numpy array

    Calls:
    -----
    percentile()

    """
    # Parameters
    depth_fraction_low = 0.4
    depth_fraction_high = 0.95
    curvature_fraction = 0.7
    high_map_value = 0.9
    precision_limit = 100

    # Take the opposite of the curvature values
    curvatures = -curvatures

#   # Normalize curvatures so that the max abs value = 1
#   curvatures = curvatures/max(np.abs(curvatures))

    # Find depth and curvature values greater than
    # the values of a fraction of the vertices
    sort_depths = np.sort(depths)
    sort_curvatures = np.sort(curvatures)
    depth_percentile_low = percentile(sort_depths, depth_fraction_low, key=lambda x:x)
    depth_percentile_high = percentile(sort_depths, depth_fraction_high, key=lambda x:x)
    curvature_percentile = percentile(sort_curvatures, curvature_fraction, key=lambda x:x)
    if verbose == 2:
        print('    depth values {0:.2f}, {1:.2f} greater than {2:.0f}%, {3:.0f}% of vertices'.
              format(depth_percentile_low, depth_percentile_high,
                     100 * depth_fraction_low, 100 * depth_fraction_high))
        print('    curvature value {0:.2f} greater than {1:.0f}% of vertices'.
              format(curvature_percentile, 100 * curvature_fraction))
        #print(sum([1. for x in depths if x < depth_percentile_low]) / sum(depths > 0))
        #print(sum([1. for x in depths if x < depth_percentile_high]) / sum(depths > 0))
        #print(sum([1. for x in curvatures if x < curvature_percentile]) / sum(depths > 0))

    # Find slope for depth and curvature values
    # (while preventing precision errors for large slope values)
    slope_factor = np.log((1. / high_map_value) - 1)
    if depth_percentile_high == depth_percentile_low:
        slope_depth = -precision_limit
        if verbose == 1:
            print('    (Warning: infinite slope depth)')
    else:
        slope_depth = slope_factor / (depth_percentile_high - depth_percentile_low)
        if slope_depth > precision_limit:
            slope_depth = precision_limit
            if verbose == 1:
                print('    (Warning: high +slope depth: ' + str(slope_depth) + ')')
        elif slope_depth < -precision_limit:
            slope_depth = -precision_limit
            if verbose == 1:
                print('    (Warning: high -slope depth: ' + str(slope_depth) + ')')

    if curvature_percentile == 0:
        slope_curvature = -precision_limit
        if verbose == 1:
            print('    (Warning: infinite slope curvature)')
    else:
        slope_curvature = slope_factor / curvature_percentile
        if slope_curvature > precision_limit:
            slope_curvature = precision_limit
            if verbose == 1:
                print('    (Warning: high +slope curvature: ' + str(slope_curvature) + ')')
        elif slope_curvature < -precision_limit:
            slope_curvature = -precision_limit
            if verbose == 1:
                print('    (Warning: high -slope curvature: ' + str(slope_curvature) + ')')

    # Map values with sigmoidal function to range [0,1]
    st_depths = 1.0 / (1.0 + np.exp(slope_depth * (depths - depth_percentile_low)))
    st_curvatures = 1.0 / (1.0 + np.exp(slope_curvature * curvatures))

    # Assign likelihood values to vertices
    likelihoods = st_depths * st_curvatures  # element-wise multiplication

    return likelihoods
