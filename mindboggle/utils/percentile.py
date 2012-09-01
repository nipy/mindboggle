#!/usr/bin/python
"""
 Compute score at percentile

 http://code.activestate.com/recipes/511478/ (r2)

 Alternative scipy implementation:
 from scipy.stats import scoreatpercentile
 depth_found = scoreatpercentile(depths, depth_threshold2)

"""

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    Parameters
    ----------
    N : list of values. Note N MUST BE already sorted
    percent : float value from 0.0 to 1.0
    key : optional key function to compute value from each element of N

    Returns
    -------
    percentile : percentile of the values

    """

    import numpy as np

    if not len(N):
        return None

    k = (len(N)-1) * percent
    f = np.floor(k)
    c = np.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)

    percentile = d0 + d1

    return percentile
