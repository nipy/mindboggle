#!/usr/bin/env python
"""
Compute functions.


Authors:
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

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

    Notes
    --------
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
        if not os.path.exists(outfile):
            raise(IOError(outfile + " not found"))
    else:
        outfile = ''

    return vector_distances, outfile

def source_to_target_distances(sourceIDs, targetIDs, points, ntargets,
                               segmentIDs=[], excludeIDs=[-1]):
    """
    Create a distance matrix between source and target points.

    Compute distances between source and target features,
    optionally within each segment.

    Example::

        Compute the minimum distance from each label boundary vertex
        (corresponding to a fundus in the DKT cortical labeling protocol)
        to all of the fundus vertices in the same fold.

    Parameters
    ----------
    sourceIDs : list of N integers (N is the number of vertices)
        source feature IDs (e.g., fundi; ignore -1)
    targetIDs : list of N integers (N is the number of vertices)
        target feature IDs (e.g., label boundaries; ignore -1)
    points : list of N lists of three floats (N is the number of vertices)
        coordinates of all vertices
    ntargets : integer
        number of total possible targets (to ensure indices are meaningful)
    segmentIDs : list of N integers (N is the number of vertices)
        segment IDs (e.g., folds; ignore -1); compute distances
        between source and target features within each segment
    excludeIDs : list of integers
        source/target/segment IDs to exclude

    Returns
    -------
    distances : numpy array
        distance value for each vertex (default -1)
    distance_matrix : numpy array [#vertices by #target features]
        distances organized by targetIDs (columns)

    """
    import numpy as np
    from mindboggle.utils.compute import point_distance

    if isinstance(points, list):
        points = np.asarray(points)

    npoints = len(points)

    distances = -1 * np.ones(npoints)
    distance_matrix = -1 * np.ones((npoints, ntargets))
    if not np.size(segmentIDs):
        segmentIDs = np.ones(len(targetIDs))

    # For each segment:
    segments = [x for x in np.unique(segmentIDs) if x not in excludeIDs]
    if segments:
        for segment in segments:
            segment_indices = [i for i,x in enumerate(segmentIDs)
                               if x == segment]

            # Find all target points in the segment:
            source_indices = [i for i,x in enumerate(sourceIDs)
                              if x not in excludeIDs
                              if i in segment_indices]
            # Find all source points in the segment:
            target_indices = [i for i,x in enumerate(targetIDs)
                              if x not in excludeIDs
                              if i in segment_indices]
            if source_indices and target_indices:

                # For each target point in the segment:
                for itarget in target_indices:

                    # Find the closest source point to the target point:
                    d, i = point_distance(points[itarget],
                                          points[source_indices])
                    distances[itarget] = d
                    distance_matrix[itarget, targetIDs[itarget] - 1] = d

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

    if np.size(W):
        X = weighted_to_repeated_values(X, W, precision)

    mad = np.median(np.abs(X - np.median(X))) / c

    return mad


def means_per_label(values, labels, exclude_labels, areas=[]):
    """
    Compute the mean value across vertices per label,
    optionally taking into account surface area per vertex.

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
    exclude_labels : list of integers
        labels to be excluded
    areas : numpy array of floats
        surface areas

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
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk
    >>> from mindboggle.utils.compute import means_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> area_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> exclude_labels = [-1]
    >>> areas = areas
    >>> #
    >>> # Example 1: compute mean curvature per label:
    >>> means, sdevs, label_list, label_areas = means_per_label(values, labels,
    >>>     exclude_labels, areas)
    >>> #
    >>> # Example 2: compute mean coordinates per label:
    >>> faces, lines, indices, points, npoints, curvs, name, input_vtk = read_vtk(values_file)
    >>> means, sdevs, label_list, label_areas = means_per_label(points, labels,
    >>>     exclude_labels, areas)

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(areas, np.ndarray):
        areas = np.asarray(areas)

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
            means.append(0)
            sdevs.append(0)
            label_areas.append(0)

    if dim > 1:
        means = [x.tolist() for x in means]
        sdevs = [x.tolist() for x in sdevs]
        label_areas = [x.tolist() for x in label_areas]

    return means, sdevs, label_list, label_areas


def sum_per_label(values, labels, exclude_labels):
    """
    Compute the sum value across vertices per label.

    Parameters
    ----------
    values : numpy array of one or more lists of integers or floats
        values to average per label
    labels : list or array of integers
        label for each value
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
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk
    >>> from mindboggle.utils.compute import sum_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> exclude_labels = [-1]
    >>> # Compute sum area per label:
    >>> sums, label_list = sum_per_label(values, labels, exclude_labels)

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)

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

def stats_per_label(values, labels, exclude_labels, weights=[], precision=1):
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
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.compute import stats_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> area_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> exclude_labels = [-1]
    >>> weights = areas
    >>> precision = 1
    >>> stats_per_label(values, labels, exclude_labels, weights, precision)

    """
    import numpy as np
    from scipy.stats import skew, kurtosis, scoreatpercentile
    from mindboggle.utils.compute import weighted_to_repeated_values, median_abs_dev

    # Make sure arguments are numpy arrays:
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(weights, np.ndarray):
        weights = np.asarray(weights)

    # Initialize all statistical lists:
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


def volume_per_label(input_file, labels=[], exclude_labels=[-1]):
    """
    Compute volume per labeled region in a nibabel-readable image.

    Parameters
    ----------
    input_file : string
        name of image file, consisting of index-labeled pixels/voxels
    labels : list of integers
        label numbers for image volumes (if empty, use unique numbers in file)
    exclude_labels : list of integers
        label IDs to exclude

    Returns
    -------
    labels_volumes : list of integer list and float list
        label numbers and volume for each labeled region (default -1)

    Examples
    --------
    >>> import os
    >>> from mindboggle.LABELS import DKTprotocol
    >>> from mindboggle.utils.compute import volume_per_label
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> input_file = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> dkt = DKTprotocol()
    >>> labels_volumes = volume_per_label(input_file, dkt.label_numbers, [-1])
    >>> print(labels_volumes)

    """
    import numpy as np
    import nibabel as nb

    # Load labeled image volumes:
    img = nb.load(input_file)
    hdr = img.get_header()
    pixdims = hdr.get_zooms()
    volume_per_voxel = np.product(pixdims)
    data = img.get_data().ravel()

    # Initialize output:
    if not labels:
        labels = np.unique(data).tolist()
    labels = [int(x) for x in labels if x not in exclude_labels]
    volumes = -1 * np.ones(len(labels))

    # Loop through labels:
    for ilabel, label in enumerate(labels):

        # Find which voxels contain the label in each volume:
        indices = np.where(data==label)[0]
        if len(indices):
            volumes[ilabel] = volume_per_voxel * len(indices)

    labels_volumes = [labels, volumes.tolist()]

    return labels_volumes


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
    >>> from mindboggle.utils.compute import compute_image_histogram
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> infile = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> compute_image_histogram(infile, nbins=100, threshold=0.1)

    """
    import numpy as np
    import nibabel as nb
    #from pylab import plot #, hist

    #-------------------------------------------------------------------------
    # Compute histogram
    #-------------------------------------------------------------------------
    # Load image
    print(infile)
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
