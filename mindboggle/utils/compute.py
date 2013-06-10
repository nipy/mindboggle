#!/usr/bin/env python
"""
Compute functions.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Compute distance
#------------------------------------------------------------------------------
def point_distance(point, points):
    """
    Compute the Euclidean distance from one point to a second (set) of points.

    Parameters
    ----------
    point : list of three floats
        coordinates for a single point
    points : list with one or more lists of three floats
        coordinates for a second point (or multiple points)

    Returns
    -------
    min_distance : float
        Euclidean distance between two points,
        or the minimum distance between a point and a set of points
    min_index : int
        index of closest of the points (zero if only one)

    Examples
    --------
    >>> from mindboggle.utils.compute import point_distance
    >>> point = [1,2,3]
    >>> points = [[10,2.0,3], [0,1.5,2]]
    >>> point_distance(point, points)
      (1.5, 1)

    """
    import numpy as np

    # If points is a single point
    if np.ndim(points) == 1:
        return np.sqrt((point[0] - points[0]) ** 2 + \
                       (point[1] - points[1]) ** 2 + \
                       (point[2] - points[2]) ** 2), 0

    # If points is a set of multiple points
    elif np.ndim(points) == 2:
        min_distance = np.Inf
        min_index = 0
        for index, point2 in enumerate(points):
            distance = np.sqrt((point[0] - point2[0]) ** 2 + \
                               (point[1] - point2[1]) ** 2 + \
                               (point[2] - point2[2]) ** 2)
            if distance < min_distance:
                min_distance = distance
                min_index = index
        return min_distance, min_index

    # Else return None
    else:
        return None, None

def vector_distance(vector1, vector2, normalize=False):
    """
    Compute the Euclidean distance between two equal-sized vectors.

    Parameters
    ----------
    vector1 : numpy array of floats
        vector of values
    vector2 : numpy array of floats
        vector of values
    normalize : Boolean
        normalize each element of the vectors?

    Returns
    -------
    distance : float
        Euclidean distance between two vectors

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.utils.compute import vector_distance
    >>> vector1 = np.array([1.,2.,3.])
    >>> vector2 = np.array([0,1,5])
    >>> vector_distance(vector1, vector2)
      0.81649658092772592

    """
    import numpy as np

    if np.size(vector1) == np.size(vector2):
    # Make sure arguments are numpy arrays
        if not isinstance(vector1, np.ndarray):
            vector1 = np.asarray(vector1)
        if not isinstance(vector2, np.ndarray):
            vector2 = np.asarray(vector2)
        if normalize:
            vector_diff = np.zeros(len(vector1))
            for i in range(len(vector1)):
                max_v1v2 = max([vector1[i], vector2[i]])
                if max_v1v2 > 0:
                    vector_diff[i] = (vector1[i] - vector2[i]) / max_v1v2
        else:
            vector_diff = vector1 - vector2
        return np.sqrt(sum((vector_diff)**2)) / np.size(vector1)
    else:
        print("Vectors have to be of equal size to compute distance.")
        return None

def pairwise_vector_distances(vectors, save_file=False, normalize=False):
    """
    Compare every pair of equal-sized vectors.

    Parameters
    ----------
    vectors : array of 1-D lists or arrays of integers or floats
    save_file : Boolean
        save file?
    normalize : Boolean
        normalize each element of the vectors?

    Returns
    -------
    vector_distances : numpy array of integers or floats
        distances between each pair of vectors
    outfile : string [optional]
        output filename for pairwise_vector_distances

    Examples
    --------
    >>> from mindboggle.utils.compute import pairwise_vector_distances
    >>> pairwise_vector_distances([[1,2,3],[0,3,5],[0,3.5,5],[1,1,1]])
        (array([[ 0.        ,  0.81649658,  0.89752747,  0.74535599],
               [ 0.        ,  0.        ,  0.16666667,  1.52752523],
               [ 0.        ,  0.        ,  0.        ,  1.60727513],
               [ 0.        ,  0.        ,  0.        ,  0.        ]]),
         '')

    """
    import os
    import numpy as np
    from mindboggle.utils.compute import vector_distance

    # Make sure argument is a numpy array
    if not isinstance(vectors, np.ndarray):
        vectors = np.array(vectors)

    # Initialize output
    vector_distances = np.zeros((len(vectors), len(vectors)))

    #---------------------------------------------------------------------------
    # Compute distance between each pair of vectors
    #---------------------------------------------------------------------------
    # Loop through every pair of vectors
    for ihist1 in range(len(vectors)):
        for ihist2 in range(len(vectors)):
            if ihist2 >= ihist1:

                # Store pairwise distances between histogram values
                d = vector_distance(1.0*vectors[ihist1],
                                            1.0*vectors[ihist2],
                                            normalize=normalize)
                vector_distances[ihist1, ihist2] = d

    if save_file:
        outfile = os.path.join(os.getcwd(), 'vector_distances.txt')
        np.savetxt(outfile, vector_distances,
                   fmt=len(vectors) * '%.4f ', delimiter='\t', newline='\n')
    else:
        outfile = ''

    return vector_distances, outfile

def weighted_to_repeated_values(X, W=[], precision=1):
    """
    Create a list of repeated values from weighted values.

    This is useful for computing weighted statistics (ex: weighted median).

    Adapted to allow for fractional weights from
        http://stackoverflow.com/questions/966896/
               code-golf-shortest-code-to-find-a-weighted-median

    Parameters
    ----------
    X : numpy array of floats or integers
        values
    W : numpy array of floats or integers
        weights
    precision : integer
        number of decimal places to consider weights

    Returns
    -------
    repeat_values : numpy array of floats or integers
        repeated values according to their weight

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.utils.compute import weighted_to_repeated_values
    >>> X = np.array([1,2,4,7,8])
    >>> W = np.array([.1,.1,.3,.2,.3])
    >>> precision = 1
    >>> weighted_to_repeated_values(X, W, precision)
        [1, 2, 4, 4, 4, 7, 7, 8, 8, 8]

    """
    import numpy as np

    # Make sure arguments have the correct type:
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(W, np.ndarray):
        W = np.array(W)
    if not isinstance(precision, int):
        precision = int(precision)

    if np.shape(W) and len(W):
        # If weights are decimals, multiply by 10 until they are whole numbers.
        # If after multiplying precision times they are not whole, round them:
        whole = True
        if any(np.mod(W,1)):
            whole = False
            for i in range(precision):
                if any(np.mod(W,1)):
                    W *= 10
                else:
                    whole = True
                    break

        if not whole:
             W = [int(np.round(x)) for x in W]

        repeat_values = sum([[x]*w for x,w in zip(X,W)],[])

    else:
        repeat_values = X

    return repeat_values

def weighted_median(X, W=[], precision=1):
    """
    Compute a weighted median.

    Parameters
    ----------
    X : numpy array of floats or integers
        values
    W : numpy array of floats or integers
        weights
    precision : integer
        number of decimal places to consider weights

    Returns
    -------
    wmedian : float
        weighted median

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.utils.compute import weighted_median
    >>> X = np.array([1,2,4,7,8])
    >>> W = np.array([.1,.1,.3,.2,.3])
    >>> precision = 1
    >>> # [1, 2, 4, 4, 4, 7, 7, 8, 8, 8]
    >>> weighted_median(X, W, precision)
        5.5

    """
    import numpy as np
    from mindboggle.utils.compute import weighted_to_repeated_values

    # Make sure arguments have the correct type:
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(W, np.ndarray):
        W = np.array(W)
    if not isinstance(precision, int):
        precision = int(precision)

    wmedian = np.median(weighted_to_repeated_values(X, W, precision))

    return wmedian

def median_abs_dev(X, W=[], precision=1, c=1.0):
    """
    Compute the (weighted) median absolute deviation.

    mad = median(abs(x - median(x))) / c

    Parameters
    ----------
    X : numpy array of floats or integers
        values
    W : numpy array of floats or integers
        weights
    precision : integer
        number of decimal places to consider weights
    c : float
        constant used as divisor for mad computation;
        c = 0.6745 is used to convert from mad to standard deviation

    Returns
    -------
    mad : float
        (weighted) median absolute deviation

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.utils.compute import median_abs_dev
    >>> X = np.array([1,2,4,7,8])
    >>> W = np.array([.1,.1,.3,.2,.3])
    >>> precision = 1
    >>> # [1, 2, 4, 4, 4, 7, 7, 8, 8, 8]
    >>> median_abs_dev(X, W, precision)
        2.0

    """
    import numpy as np
    from mindboggle.utils.compute import weighted_to_repeated_values

    # Make sure arguments have the correct type:
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(W, np.ndarray):
        W = np.array(W)
    if not isinstance(precision, int):
        precision = int(precision)

    if np.shape(W) and len(W):
        X = weighted_to_repeated_values(X, W, precision)

    mad = np.median(np.abs(X - np.median(X))) / c

    return mad

