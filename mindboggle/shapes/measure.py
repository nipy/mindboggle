#!/usr/bin/env python
"""
Shape calculations.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net

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
    >>> from mindboggle.shapes.measure import point_distance
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
    >>> from mindboggle.shapes.measure import vector_distance
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
    >>> from mindboggle.shapes.measure import pairwise_vector_distances
    >>> pairwise_vector_distances([[1,2,3],[0,3,5],[0,3.5,5],[1,1,1]])
        (array([[ 0.        ,  0.81649658,  0.89752747,  0.74535599],
               [ 0.        ,  0.        ,  0.16666667,  1.52752523],
               [ 0.        ,  0.        ,  0.        ,  1.60727513],
               [ 0.        ,  0.        ,  0.        ,  0.        ]]),
         '')

    """
    import os
    import numpy as np
    from mindboggle.shapes.measure import vector_distance

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

def area(command, surface_file):
    """
    Measure area of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    Parameters
    ----------
    command : Voronoi-based surface area C++ executable command
    surface_file : ``vtk file``

    """
    import os
    from nipype.interfaces.base import CommandLine

    area_file = os.path.join(os.getcwd(),
                os.path.splitext(os.path.basename(surface_file))[0] + '.area.vtk')
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, area_file])
    cli.cmdline
    cli.run()

    return area_file

def depth(command, surface_file):
    """
    Measure "travel depth" of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    Parameters
    ----------
    command : travel depth C++ executable command
    surface_file : ``vtk file``

    """
    import os
    from nipype.interfaces.base import CommandLine

    depth_file = os.path.join(os.getcwd(),
                 os.path.splitext(os.path.basename(surface_file))[0] + '.depth.vtk')
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, depth_file])
    cli.cmdline
    cli.run()

    return depth_file

def curvature(command, surface_file):
    """
    Measure curvature values of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    command : curvature C++ executable command
    surface_file : ``vtk file``

    """
    import os
    from nipype.interfaces.base import CommandLine

    stem = os.path.join(os.getcwd(),
                        os.path.splitext(os.path.basename(surface_file))[0])
    mean_curvature_file = stem + '.curv.avg.vtk'
    gauss_curvature_file = stem + '.curv.gauss.vtk'
    max_curvature_file = stem + '.curv.max.vtk'
    min_curvature_file = stem + '.curv.min.vtk'
    min_curvature_vector_file = stem + '.curv.min.dir.txt'
    args = ['-g', gauss_curvature_file,
            '-x', max_curvature_file,
            '-i', min_curvature_file,
            '-d', min_curvature_vector_file,
            surface_file, mean_curvature_file]
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    return mean_curvature_file, gauss_curvature_file,\
           max_curvature_file, min_curvature_file, min_curvature_vector_file

def mean_value_per_label(values, areas, labels, exclude_labels):
    """
    Compute the mean value across vertices per label,
    taking into account surface area per vertex.

    average value = sum(a_i * v_i) / total_surface_area,
    where *a_i* and *v_i* are the area and value for each vertex *i*.

    Parameters
    ----------
    values : numpy array of integer or float values
    areas : numpy array of surface areas
    labels : list or array of integer labels (same length as values)
    exclude_labels : list of integer labels to be excluded

    Returns
    -------
    mean_values : list of floats
        mean values
    norm_mean_values : list of floats
        mean values normalized by vertex area
    surface_areas : list of floats
        surface area for each labeled set of vertices
    label_list : list of integers
        unique label numbers

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.shapes.measure import mean_value_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'arno', 'measures', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> depths, name = read_scalars(depth_file, True, True)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> exclude_labels = [-1,0]
    >>> mean_values, norm_mean_values, surface_areas, label_list = mean_value_per_label(depths,
    >>>     areas, labels, exclude_labels)

    """
    import numpy as np

    def avg_by_area(values_label, areas_label):
        return sum(areas_label * values_label) / sum(areas_label)

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(areas, np.ndarray):
        areas = np.asarray(areas)

    label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]
    mean_values = []
    norm_mean_values = []
    surface_areas = []

    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            mean_value = np.mean(values[I])
            norm_mean_value = avg_by_area(values[I], areas[I])
            mean_values.append(mean_value)
            norm_mean_values.append(norm_mean_value)
            surface_area = sum(areas[I])
            surface_areas.append(surface_area)
        else:
            mean_values.append(0)
            norm_mean_values.append(0)
            surface_areas.append(0)

    return mean_values, norm_mean_values, surface_areas, label_list

def volume_per_label(labels, input_file):
    """
    Compute volume per labeled region in a nibabel-readable (e.g., nifti) image.

    Parameters
    ----------
    labels : list of integers
        label numbers for image volumes
    input_file : string
        name of image file, consisting of index-labeled pixels/voxels

    Returns
    -------
    volumes : list of floats
        volume for each labeled region
    labels : list of integers
        label numbers for image volumes

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_file import read_columns
    >>> from mindboggle.shapes.measure import volume_per_label
    >>> path = os.path.join(os.environ['MINDBOGGLE_DATA'])
    >>> input_file = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> labels_file = os.path.join(path, 'info', 'labels.volume.DKT25.txt')
    >>> labels = read_columns(labels_file, 1)[0]
    >>> volumes = volume_per_label(labels, input_file)
    >>> print(volumes)

    """
    import numpy as np
    import nibabel as nb

    # Load labeled image volumes
    data = nb.load(input_file).get_data().ravel()

    # Initialize output
    volumes = np.zeros((len(labels), 1))

    # Loop through labels
    for ilabel, label in enumerate(labels):
        label = int(label)
        volumes[ilabel, 0] = label

        # Find which voxels contain the label in each volume
        indices = np.where(data==label)[0]
        volumes[ilabel] = len(indices)

    return volumes, labels

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    http://code.activestate.com/recipes/511478/ (r2)

    Alternative scipy implementation:
    from scipy.stats import scoreatpercentile
    depth_found = scoreatpercentile(depths, depth_threshold2)

    Parameters
    ----------
    N : list of values. Note N MUST BE already sorted
    percent : float value from 0.0 to 1.0
    key : optional key function to compute value from each element of N

    Returns
    -------
    percentile : percentile of the values

    Examples
    --------
    >>> from mindboggle.shapes.measure import percentile
    >>> N = [2,3,4,8,9,10]
    >>> percent = 0.5
    >>> percentile(N, percent)
      6.0

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
