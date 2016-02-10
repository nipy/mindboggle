#!/usr/bin/env python
"""
Compute functions.


Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


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
    >>> from mindboggle.guts.compute import point_distance
    >>> point = [1,2,3]
    >>> points = [[10,2.0,3], [0,1.5,2]]
    >>> point_distance(point, points)
    (1.5, 1)

    Notes
    -----
    Future plan is to use scipy.spatial.distance.cdist to compute distances
    scipy.spatial.distance.cdist is available in scipy v0.12 or later

    """
    import numpy as np

    # If points is a single point
    if np.ndim(points) == 1:
       #return np.linalg.norm(np.array(point) - np.array(points))
       return np.sqrt((point[0] - points[0]) ** 2 + \
                       (point[1] - points[1]) ** 2 + \
                       (point[2] - points[2]) ** 2), 0

    # If points is a set of multiple points
    elif np.ndim(points) == 2:
        min_distance = np.Inf
        min_index = 0
        point = np.array(point)
        for index, point2 in enumerate(points):
            distance = np.sqrt((point[0] - point2[0]) ** 2 + \
                               (point[1] - point2[1]) ** 2 + \
                               (point[2] - point2[2]) ** 2)

            #distance = np.linalg.norm(point - np.array(point2))

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
    >>> from mindboggle.guts.compute import vector_distance
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
    >>> from mindboggle.guts.compute import pairwise_vector_distances
    >>> vectors = [[1,2,3],[0,3,5],[0,3.5,5],[1,1,1]]
    >>> save_file = False
    >>> normalize = False
    >>> vector_distances, outfile = pairwise_vector_distances(vectors,
    ...     save_file, normalize)
    >>> vector_distances
    array([[ 0.        ,  0.81649658,  0.89752747,  0.74535599],
           [ 0.        ,  0.        ,  0.16666667,  1.52752523],
           [ 0.        ,  0.        ,  0.        ,  1.60727513],
           [ 0.        ,  0.        ,  0.        ,  0.        ]])

    """
    import os
    import numpy as np
    from mindboggle.guts.compute import vector_distance

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
        if not os.path.exists(outfile):
            raise IOError(outfile + " not found")
    else:
        outfile = ''

    return vector_distances, outfile


def source_to_target_distances(sourceIDs, targetIDs, points,
                               segmentIDs=[], excludeIDs=[-1]):
    """
    Create a Euclidean distance matrix between source and target points.

    Compute the Euclidean distance from each source point to
    its nearest target point, optionally within each segment.

    Example::

        Compute fundus-to-feature distances, the minimum distance
        from each label boundary vertex (corresponding to a fundus
        in the DKT cortical labeling protocol) to all of the
        feature vertices in the same fold.

    Parameters
    ----------
    sourceIDs : list of N integers (N is the number of vertices)
        source IDs, where any ID not in excludeIDs is a source point
    targetIDs : list of N integers (N is the number of vertices)
        target IDs, where any ID not in excludeIDs is a target point
    points : list of N lists of three floats (N is the number of vertices)
        coordinates of all vertices
    segmentIDs : list of N integers (N is the number of vertices)
        segment IDs, where each ID not in excludeIDs is considered a
        different segment (unlike above, where value in sourceIDs or
        targetIDs doesn't matter, so long as its not in excludeIDs);
        source/target distances are computed within each segment
    excludeIDs : list of integers
        IDs to exclude

    Returns
    -------
    distances : numpy array
        distance value for each vertex (default -1)
    distance_matrix : numpy array [#points by maximum segment ID + 1]
        distances organized by segments (columns)

    """
    import numpy as np
    from mindboggle.guts.compute import point_distance

    if isinstance(points, list):
        points = np.asarray(points)
    npoints = len(points)

    # Extract unique segment IDs (or use all points as a single segment):
    if np.size(segmentIDs):
        segments = [x for x in np.unique(segmentIDs) if x not in excludeIDs]
    else:
        segmentIDs = np.zeros(npoints)
        segments = [0]
    nsegments = max(segments) + 1

    # Initialize outputs:
    distances = -1 * np.ones(npoints)
    distance_matrix = -1 * np.ones((npoints, nsegments))

    # For each segment:
    for segment in segments:
        segment_indices = [i for i,x in enumerate(segmentIDs)
                           if x == segment]

        # Find all source points in the segment:
        source_indices = [i for i,x in enumerate(sourceIDs)
                          if x not in excludeIDs
                          if i in segment_indices]
        # Find all target points in the segment:
        target_indices = [i for i,x in enumerate(targetIDs)
                          if x not in excludeIDs
                          if i in segment_indices]

        if source_indices and target_indices:

            # For each source point in the segment:
            for isource in source_indices:

                # Find the closest target point:
                d, i = point_distance(points[isource],
                                      points[target_indices])
                distances[isource] = d
                distance_matrix[isource, segment] = d

    return distances, distance_matrix


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
    >>> from mindboggle.guts.compute import weighted_to_repeated_values
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

    if np.size(W):
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
    >>> from mindboggle.guts.compute import weighted_median
    >>> X = np.array([1,2,4,7,8])
    >>> W = np.array([.1,.1,.3,.2,.3])
    >>> precision = 1
    >>> # [1, 2, 4, 4, 4, 7, 7, 8, 8, 8]
    >>> weighted_median(X, W, precision)
    5.5

    """
    import numpy as np
    from mindboggle.guts.compute import weighted_to_repeated_values

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
    >>> from mindboggle.guts.compute import median_abs_dev
    >>> X = np.array([1,2,4,7,8])
    >>> W = np.array([.1,.1,.3,.2,.3])
    >>> precision = 1
    >>> # [1, 2, 4, 4, 4, 7, 7, 8, 8, 8]
    >>> median_abs_dev(X, W, precision)
    2.0

    """
    import numpy as np
    from mindboggle.guts.compute import weighted_to_repeated_values

    # Make sure arguments have the correct type:
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(W, np.ndarray):
        W = np.array(W)
    if not isinstance(precision, int):
        precision = int(precision)

    if np.size(W):
        X = weighted_to_repeated_values(X, W, precision)

    mad = np.median(np.abs(X - np.median(X))) / c

    return mad


def means_per_label(values, labels, include_labels=[], exclude_labels=[], areas=[]):
    """
    Compute the mean value across vertices per label,
    optionally taking into account surface area per vertex (UNTESTED).

    Formula:
    average value = sum(a_i * v_i) / total_surface_area,
    where *a_i* and *v_i* are the area and value for each vertex *i*.

    Note ::
        This function is different than stats_per_label() in two ways:
            1. It only computes the (weighted) mean and sdev.
            2. It can accept 2-D arrays (such as [x,y,z] coordinates).

    Parameters
    ----------
    values : numpy array of one or more lists of integers or floats
        values to average per label
    labels : list or array of integers
        label for each value
    include_labels : list of integers
        labels to include
    exclude_labels : list of integers
        labels to be excluded
    areas : numpy array of floats
        surface areas (if provided, used to normalize means and sdevs)

    Returns
    -------
    means : list of floats
        mean(s) for each label
    sdevs : list of floats
        standard deviation(s) for each label
    label_list : list of integers
        unique label numbers
    label_areas : list of floats (if normalize_by_area)
        surface area for each labeled set of vertices

    Examples
    --------
    >>> from mindboggle.mio.vtks import read_scalars, read_vtk
    >>> from mindboggle.guts.compute import means_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> values_file = fetch_data(urls['left_mean_curvature'])
    >>> labels_file = fetch_data(urls['left_freesurfer_labels'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> values, name = read_scalars(values_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> include_labels = []
    >>> exclude_labels = [-1]
    >>> # Compute mean curvature per label normalized by area:
    >>> means, sdevs, label_list, label_areas = means_per_label(values, labels,
    ...     include_labels, exclude_labels, areas)
    >>> means[0:3]
    [-1.1793044671582817, -1.2140542483504968, -2.493175052379216]
    >>> sdevs[0:3]
    [2.4382673224633384, 2.3385730827461524, 2.0185030879621064]
    >>> # Compute mean curvature per label:
    >>> areas = []
    >>> means, sdevs, label_list, label_areas = means_per_label(values, labels,
    ...     include_labels, exclude_labels, areas)
    >>> means[0:3]
    [-0.99076805512931965, -0.30049554123055167, -1.5934202292058071]
    >>> sdevs[0:3]
    [2.3486044266667516, 2.4022970537134016, 2.3252993828800501]

    >>> # FIX: compute mean coordinates per label:
    >>> #points, indices, lines, faces, labels, scalar_names, npoints, input_vtk = read_vtk(values_file)
    >>> #means, sdevs, label_list, label_areas = means_per_label(points, labels,
    >>> #    include_labels, exclude_labels, areas)
    >>> #means[0:3]
    >>> #sdevs[0:3]

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(areas, np.ndarray):
        areas = np.asarray(areas)

    if include_labels:
        label_list = include_labels
    else:
        label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]
    means = []
    sdevs = []
    label_areas = []
    if values.ndim > 1:
        dim = np.shape(values)[1]
    else:
        dim = 1

    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            X = values[I]
            if np.size(areas):
                W = areas[I]
                sumW = sum(W)
                label_areas.append(sumW)
                if sumW > 0:
                    if dim > 1:
                        W = np.transpose(np.tile(W, (dim,1)))
                    means.append(np.sum(W * X, axis=0) / sumW)
                    Xdiff = X - np.mean(X, axis=0)
                    sdevs.append(np.sqrt(np.sum(W * Xdiff**2, axis=0) / sumW))
                else:
                    if dim > 1:
                        means.append(np.mean(X, axis=0))
                        sdevs.append(np.std(X, axis=0))
                    else:
                        means.append(np.mean(X))
                        sdevs.append(np.std(X))
            else:
                if dim > 1:
                    means.append(np.mean(X, axis=0))
                    sdevs.append(np.std(X, axis=0))
                else:
                    means.append(np.mean(X))
                    sdevs.append(np.std(X))
        else:
            means.append(np.zeros(dim))
            sdevs.append(np.zeros(dim))
            label_areas.append(np.zeros(dim))

    if dim > 1:
        means = [x.tolist() for x in means]
        sdevs = [x.tolist() for x in sdevs]
        label_areas = [x.tolist() for x in label_areas]

    return means, sdevs, label_list, label_areas


def sum_per_label(values, labels, include_labels=[], exclude_labels=[]):
    """
    Compute the sum value across vertices per label.

    Parameters
    ----------
    values : numpy array of one or more lists of integers or floats
        values to average per label
    labels : list or array of integers
        label for each value
    include_labels : list of integers
        labels to include
    exclude_labels : list of integers
        labels to be excluded

    Returns
    -------
    sums : list of floats
        sum for each label
    label_list : list of integers
        unique label numbers

    Examples
    --------
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.compute import sum_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> values_file = fetch_data(urls['left_mean_curvature'])
    >>> labels_file = fetch_data(urls['left_freesurfer_labels'])
    >>> values, name = read_scalars(values_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> include_labels = []
    >>> exclude_labels = [-1]
    >>> # Compute sum area per label:
    >>> sums, label_list = sum_per_label(values, labels, include_labels,
    ...                                  exclude_labels)
    >>> sums[0:3]
    [-8228.3286978489996, -424.90069530000005, -1865.8950884000001]

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)

    if include_labels:
        label_list = include_labels
    else:
        label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]
    sums = []
    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            X = values[I]
            sums.append(np.sum(X))
        else:
            sums.append(0)

    return sums, label_list


def stats_per_label(values, labels, include_labels=[], exclude_labels=[],
                    weights=[], precision=1):
    """
    Compute various statistical measures across vertices per label,
    optionally using weights (such as surface area per vertex).

    Example (area-weighted mean):
    average value = sum(a_i * v_i) / total_surface_area,
    where *a_i* and *v_i* are the area and value for each vertex *i*.

    Note ::
        This function is different than means_per_label() in two ways:
            1. It computes more than simply the (weighted) mean and sdev.
            2. It only accepts 1-D arrays of values.

    Reference
    ---------
    Weighted skewness and kurtosis unbiased by sample size
    Lorenzo Rimoldini, arXiv:1304.6564 (2013)
    http://arxiv.org/abs/1304.6564

    Parameters
    ----------
    values : numpy array of individual or lists of integers or floats
        values for all vertices
    labels : list or array of integers
        label for each value
    include_labels : list of integers
        labels to include
    exclude_labels : list of integers
        labels to be excluded
    weights : numpy array of floats
        weights to compute weighted statistical measures
    precision : integer
        number of decimal places to consider weights

    Returns
    -------
    medians : list of floats
        median for each label
    mads : list of floats
        median absolute deviation for each label
    means : list of floats
        mean for each label
    sdevs : list of floats
        standard deviation for each label
    skews : list of floats
        skew for each label
    kurts : list of floats
        kurtosis value for each label
    lower_quarts : list of floats
        lower quartile for each label
    upper_quarts : list of floats
        upper quartile for each label
    label_list : list of integers
        list of unique labels

    Examples
    --------
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.compute import stats_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> values_file = fetch_data(urls['left_mean_curvature'])
    >>> labels_file = fetch_data(urls['left_freesurfer_labels'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> values, name = read_scalars(values_file, True, True)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> include_labels = []
    >>> exclude_labels = [-1]
    >>> weights = areas
    >>> precision = 1
    >>> medians, mads, means, sdevs, skews, kurts, lower_quarts, upper_quarts, label_list = stats_per_label(values,
    ...     labels, include_labels, exclude_labels, weights, precision)
    >>> medians[0:3]
    [-1.13602, -1.2296100000000001, -2.4966499999999998]
    >>> mads[0:3]
    [1.170255, 1.5045030000000001, 1.28234]
    >>> means[0:3]
    [-1.1793044671582793, -1.2140542483504964, -2.4931750523792142]
    >>> kurts[0:3]
    [2.3411784448846227, -0.39690083155055289, -0.55786762009358526]

    """
    import numpy as np
    from scipy.stats import skew, kurtosis, scoreatpercentile
    from mindboggle.guts.compute import weighted_to_repeated_values, median_abs_dev

    # Make sure arguments are numpy arrays:
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(weights, np.ndarray):
        weights = np.asarray(weights)

    # Initialize all statistical lists:
    if include_labels:
        label_list = include_labels
    else:
        label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]
    medians = []
    mads = []
    means = []
    sdevs = []
    skews = []
    kurts = []
    lower_quarts = []
    upper_quarts = []

    # Extract all vertex indices for each label:
    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            # Get the vertex values:
            X = values[I]
            if len([x for x in X if x != 0]):
                # If there are as many weights as values, apply the weights to the values:
                if np.size(weights) == np.size(values):
                    W = weights[I]
                    sumW = np.sum(W)
                    # If the sum of the weights and the standard deviation is non-zero,
                    # compute all statistics of the weighted values:
                    if sumW > 0:
                        Xdiff = X - np.mean(X)
                        Xstd = np.sqrt(np.sum(W * Xdiff**2) / sumW)
                        means.append(np.sum(W * X) / sumW)
                        sdevs.append(Xstd)
                        if Xstd > 0:
                            skews.append((np.sum(W * Xdiff**3) / sumW) / Xstd**3)
                            kurts.append((np.sum(W * Xdiff**4) / sumW) / Xstd**4 - 3)
                        else:
                            skews.append(skew(X))
                            kurts.append(kurtosis(X))
                        X = weighted_to_repeated_values(X, W, precision)
                    # If the sum of the weights equals zero, simply compute the statistics:
                    else:
                        means.append(np.mean(X))
                        sdevs.append(np.std(X))
                        skews.append(skew(X))
                        kurts.append(kurtosis(X))
                # If there are no (or not enough) weights, simply compute the statistics:
                else:
                    means.append(np.mean(X))
                    sdevs.append(np.std(X))
                    skews.append(skew(X))
                    kurts.append(kurtosis(X))
                # Compute median, median absolute deviation, and lower and upper quartiles:
                if np.size(X):
                    medians.append(np.median(X))
                    mads.append(median_abs_dev(X))
                    lower_quarts.append(scoreatpercentile(X, 25))
                    upper_quarts.append(scoreatpercentile(X, 75))
                # If the weights are all smaller than the precision, then X will disappear,
                # so set the above statistics (in the 'if' block) to zero:
                else:
                    medians.append(0)
                    mads.append(0)
                    lower_quarts.append(0)
                    upper_quarts.append(0)
            # If all values are equal to zero, set all statistics to zero:
            else:
                medians.append(0)
                mads.append(0)
                means.append(0)
                sdevs.append(0)
                skews.append(0)
                kurts.append(0)
                lower_quarts.append(0)
                upper_quarts.append(0)
        # If there are no vertices for the label, set all statistics to zero:
        else:
            medians.append(0)
            mads.append(0)
            means.append(0)
            sdevs.append(0)
            skews.append(0)
            kurts.append(0)
            lower_quarts.append(0)
            upper_quarts.append(0)

    return medians, mads, means, sdevs, skews, kurts, \
           lower_quarts, upper_quarts, label_list


def count_per_label(labels, include_labels=[], exclude_labels=[]):
    """
    Compute the number of times each label occurs.

    Parameters
    ----------
    labels : numpy 1-D array of integers
        labels (e.g., one label per vertex of a mesh)
    include_labels : list of integers
        labels to include
        (if empty, use unique numbers in image volume file)
    exclude_labels : list of integers
        labels to be excluded

    Returns
    -------
    unique_labels : list of integers
        unique label numbers
    counts : list of floats
        number of times each label occurs

    Examples
    --------
    >>> from mindboggle.guts.compute import count_per_label
    >>> labels = [8,8,8,8,8,10,11,12,10,10,11,11,11,12,12,12,12,13]
    >>> include_labels = [9,10,11,12]
    >>> exclude_labels = [13]
    >>> unique_labels, counts = count_per_label(labels, include_labels,
    ...                                         exclude_labels)
    >>> unique_labels
    [10, 11, 12]
    >>> counts
    [3, 4, 5]

    Skip the following when testing on https://travis-ci.org/nipy/mindboggle
    because the travis.yml file doesn't install nibabel (requires sudo):

    >>> import nibabel as nb  # doctest: +SKIP
    >>> from mindboggle.mio.vtks import read_scalars  # doctest: +SKIP
    >>> from mindboggle.mio.labels import DKTprotocol  # doctest: +SKIP
    >>> from mindboggle.guts.compute import count_per_label  # doctest: +SKIP
    >>> from mindboggle.mio.fetch_data import prep_tests  # doctest: +SKIP
    >>> urls, fetch_data = prep_tests()  # doctest: +SKIP
    >>> labels_file = fetch_data(urls['freesurfer_labels'], 'test.nii.gz')  # doctest: +SKIP
    >>> img = nb.load(labels_file)  # doctest: +SKIP
    >>> hdr = img.get_header()  # doctest: +SKIP
    >>> labels = img.get_data().ravel()  # doctest: +SKIP
    >>> dkt = DKTprotocol()  # doctest: +SKIP
    >>> include_labels = dkt.label_numbers  # doctest: +SKIP
    >>> exclude_labels = []  # doctest: +SKIP
    >>> unique_labels, counts = count_per_label(labels,
    ...     include_labels, exclude_labels)  # doctest: +SKIP
    >>> counts[0:5]  # doctest: +SKIP
    [972, 2414, 2193, 8329, 2941]

    """
    import numpy as np

    # Make sure labels is a numpy array:
    if isinstance(labels, list):
        labels = np.array(labels)
    elif isinstance(labels, ndarray):
        pass
    else:
        raise IOError("labels should be a numpy array.")

    # Unique list of labels:
    if include_labels:
        label_list = include_labels
    else:
        label_list = np.unique(labels).tolist()
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]

    # Loop through labels:
    unique_labels = []
    counts = []
    for ilabel, label in enumerate(label_list):

        # Find which voxels contain the label in each volume:
        indices = np.where(labels == label)[0]
        if len(indices):
            count = len(indices)
            unique_labels.append(label)
            counts.append(count)

    return unique_labels, counts


def compute_overlaps(targets, list1, list2, output_file='', save_output=True,
                     verbose=False):
    """
    Compute overlap for each target between two lists of numbers.

    Parameters
    ----------
    targets : list of integers
        targets to find in lists
    list1 : 1-D numpy array (or list) of numbers
    list2 : 1-D numpy array (or list) of numbers
    output_file : string
        (optional) output file name
    save_output : Boolean
        save output file?
    verbose : Boolean
        print statements?

    Returns
    -------
    dice_overlaps : numpy array
        Dice overlap values
    jacc_overlaps : numpy array
        Jaccard overlap values
    output_file : string
        output text file name with overlap values

    Examples
    --------
    >>> from mindboggle.guts.compute import compute_overlaps
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> list1 = [1021, 1021, 1021, 1021, 1021, 1010, 1010, 1010, 1010, 1010]
    >>> list2 = [1003, 1021, 1021, 1021, 1021, 1021, 1003, 1003, 1003, 1003]
    >>> dkt = DKTprotocol()
    >>> targets = dkt.cerebrum_cortex_DKT31_numbers
    >>> output_file = ''
    >>> save_output = True
    >>> verbose = False
    >>> dice_overlaps, jacc_overlaps, output_file = compute_overlaps(targets,
    ...     list1, list2, output_file, save_output, verbose)
    >>> dice_overlaps[18]
    0.80000000000000004
    >>> jacc_overlaps[18]
    0.66666666666666663

    """
    import os
    import numpy as np
    import pandas as pd

    if isinstance(list1, list):
        list1 = np.array(list1)
    if isinstance(list2, list):
        list2 = np.array(list2)

    if np.size(list1) != np.size(list2):
        raise IOError("Files are different sizes")

    # Initialize output:
    dice_overlaps = np.zeros(len(targets))
    jacc_overlaps = np.zeros(len(targets))
    if save_output and not output_file:
        output_file = os.path.join(os.getcwd(), 'ID_dice_jaccard.csv')

    # Loop through targets:
    for itarget, target in enumerate(targets):

        # Find which indices contain the target:
        list1_indices = np.where(list1 == target)[0]
        list2_indices = np.where(list2 == target)[0]
        len1 = len(list1_indices)
        len2 = len(list2_indices)

        # Compute their intersection and union:
        len_intersection = len(np.intersect1d(list2_indices, list1_indices))
        len_union = len(np.union1d(list2_indices, list1_indices))

        # If there is at least one target in each list:
        if len2 * len1 > 0:

            # Compute Dice and Jaccard coefficients:
            dice = np.float(2.0 * len_intersection) / (len2 + len1)
            jacc = np.float(len_intersection) / len_union
            dice_overlaps[itarget] = dice
            jacc_overlaps[itarget] = jacc
            if verbose:
                print('target: {0}, dice: {1:.2f}, jacc: {2:.2f}'.format(
                      target, dice, jacc))

    # Save output:
    if save_output:
        #np.savetxt(output_file, overlaps, fmt='%d %.4f %.4f',
        #           delimiter='\t', newline='\n')
        df1 = pd.DataFrame({'ID': targets})
        df2 = pd.DataFrame({'Dice overlap': dice_overlaps})
        df3 = pd.DataFrame({'Jaccard overlap': jacc_overlaps})
        df = pd.concat([df1, df2, df3], axis=1)
        df.to_csv(output_file, index=False)

    return dice_overlaps, jacc_overlaps, output_file


def compute_image_histogram(infile, nbins=100, threshold=0.0):
    """
    Compute histogram values from nibabel-readable image.

    Parameters
    ----------
    infile : string
        input nibabel-readable image file name
    nbins : integer
        number of bins
    threshold : float
        remove values lower than threshold

    Returns
    -------
    histogram_values : numpy array
        histogram bin values

    Examples
    --------
    >>> import os
    >>> from mindboggle.guts.compute import compute_image_histogram
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> labels_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(labels_file, labels_file + '.nii.gz')
    >>> labels_file = labels_file + '.nii.gz'
    >>> nbins = 100
    >>> threshold = 0.5
    >>> histogram_values = compute_image_histogram(labels_file, nbins,
    ...                                            threshold)
    >>> histogram_values[0:3]
    array([102865, 119610,      0])

    """
    import numpy as np
    import nibabel as nb
    #from pylab import plot #, hist

    #-------------------------------------------------------------------------
    # Compute histogram
    #-------------------------------------------------------------------------
    # Load image
    data = nb.load(infile).get_data().ravel()

    # Threshold image
    if threshold > 0:
        data = data / max(data)
        data = data[data >= threshold]

    # Compute histogram
    histogram_values, bin_edges = np.histogram(data, bins=nbins)

    # plot(range(len(histogram_values)), histogram_values, '-')
    ##a,b,c = hist(data, bins=nbins)

    return histogram_values


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()
