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
    precision_limit = 100
    #high_map_value = 0.9  # inserted below to reduce computation of the slope_gain_factor

    #=======================
    # Reset curvature values
    #=======================
    # Take the opposite of the curvature values
    curvatures = -curvatures

    # Normalize curvatures to the interval [0,1]
    curvatures -= min(curvatures)
    curvatures /= max(curvatures)

    #=============================================
    # Find depth and curvature values greater than
    # the values of a fraction of the vertices
    #=============================================
    sort_depths = np.sort(depths)
    sort_curvatures = np.sort(curvatures)
    depth_percentile_low = percentile(sort_depths, depth_fraction_low, key=lambda x:x)
    depth_percentile_high = percentile(sort_depths, depth_fraction_high, key=lambda x:x)
    curvature_percentile = percentile(sort_curvatures, curvature_fraction, key=lambda x:x)
    if verbose == 2:
        print('    depth values {:.2f}, {:.2f} greater than {:.0f}%, {:.0f}% of vertices'.
              format(depth_percentile_low, depth_percentile_high,
                     100 * depth_fraction_low, 100 * depth_fraction_high))
        print('    curvature value {:.2f} greater than {:.0f}% of vertices'.
              format(curvature_percentile, 100 * curvature_fraction))
        #print(sum([1. for x in depths if x < depth_percentile_low]) / sum(depths > 0))
        #print(sum([1. for x in depths if x < depth_percentile_high]) / sum(depths > 0))
        #print(sum([1. for x in curvatures if x < curvature_percentile]) / sum(depths > 0))

    #==========================================
    # Find slope for depth and curvature values
    #==========================================
    # Slope factor for "gain" or "sharpness" of the sigmoidal function below
    slope_factor = -2.1972245773362191  #np.log((1. / high_map_value) - 1)

    # If the slope denominator would not equal zero, compute slope
    if depth_percentile_high != depth_percentile_low:

        #----------------------------
        # Find slope for depth values
        #----------------------------
        slope_gain_depth = slope_factor / (depth_percentile_high - depth_percentile_low)

        # Prevent precision errors for large slope values
        if slope_gain_depth > precision_limit:
            slope_gain_depth = precision_limit
            if verbose == 1:
                print('    (Warning: high +slope depth)')
        elif slope_gain_depth < -precision_limit:
            slope_gain_depth = -precision_limit
            if verbose == 1:
                print('    (Warning: high -slope depth)')
    # If the denominator is equal to zero, set slope
    else:
        slope_gain_depth = -precision_limit
        if verbose == 1:
            print('    (Warning: infinite slope depth)')

    # If the slope denominator would not equal zero, compute slope
    if curvature_percentile != 0:

        #--------------------------------
        # Find slope for curvature values
        #--------------------------------
        slope_gain_curvature = slope_factor / curvature_percentile

        # Prevent precision errors for large slope values
        if slope_gain_curvature > precision_limit:
            slope_gain_curvature = precision_limit
            if verbose == 1:
                print('    (Warning: high +slope curvature)')
        elif slope_gain_curvature < -precision_limit:
            slope_gain_curvature = -precision_limit
            if verbose == 1:
                print('    (Warning: high -slope curvature)')
    # If the denominator is equal to zero, set slope
    else:
        slope_gain_curvature = -precision_limit
        if verbose == 1:
            print('    (Warning: infinite slope curvature)')

    #==========================
    # Compute likelihood values
    #==========================
    # Map values with sigmoidal function to range [0,1]
    # Y(t) = 1/(1 + e(k*(X - thr)), where k is the gain factor, and thr is the slider term
    sigmoidal_depths = 1 / (1.0 + np.exp(slope_gain_depth * (depths - depth_percentile_low)))
    sigmoidal_curvatures = 1 / (1.0 + np.exp(slope_gain_curvature * curvatures))

    # Assign likelihood values to vertices
    likelihoods = sigmoidal_depths * sigmoidal_curvatures  # element-wise multiplication

    return likelihoods
