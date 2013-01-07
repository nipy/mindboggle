#!/usr/bin/python

"""
Shape calculations.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import os
#import numpy as np
#from nipype.interfaces.base import CommandLine

##############################################################################
#   Surface calculations
##############################################################################

#------------------------------------------------------------------------------
# Compute distance
#------------------------------------------------------------------------------
def compute_point_distance(point, points):
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
    >>> from mindboggle.measure.measure_functions import compute_point_distance
    >>> point = [1,2,3]
    >>> points = [[10,2.0,3], [0,1.5,2]]
    >>> compute_point_distance(point, points)
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

def compute_vector_distance(vector1, vector2, normalize=False):
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
    >>> from mindboggle.measure.measure_functions import compute_vector_distance
    >>> vector1 = np.array([1.,2.,3.])
    >>> vector2 = np.array([0,1,5])
    >>> compute_vector_distance(vector1, vector2)
      0.81649658092772592

    """
    import numpy as np

    if np.size(vector1) == np.size(vector2):
        if type(vector1) != np.ndarray:
            vector1 = np.asarray(vector1)
        if type(vector2) != np.ndarray:
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
    vectors : list or array of 1-D lists or arrays of integers or floats
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
    >>> import numpy as np
    >>> from mindboggle.measure.measure_functions import pairwise_vector_distances
    >>> pairwise_vector_distances([[1,2,3],[0,3,5],[0,3.5,5],[1,1,1]])
        (array([[ 0.        ,  0.81649658,  0.89752747,  0.74535599],
               [ 0.        ,  0.        ,  0.16666667,  1.52752523],
               [ 0.        ,  0.        ,  0.        ,  1.60727513],
               [ 0.        ,  0.        ,  0.        ,  0.        ]]),
         '')

    """
    import os
    import numpy as np
    from mindboggle.measure.measure_functions import compute_vector_distance

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
                d = compute_vector_distance(1.0*vectors[ihist1],
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

def compute_area(command, surface_file):
    """
    Measure Joachim Giard's "travel area" for a surface mesh.

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

def compute_depth(command, surface_file):
    """
    Measure Joachim Giard's "travel depth" for a surface mesh.

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

def compute_curvature(command, surface_file):
    """
    Measure curvatures for a surface mesh.

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

def mean_value_per_label(values, areas, labels, exclude_values):
    """
    Compute the mean value across vertices per label,
    taking into account surface area per vertex.

    average value = sum(a_i * v_i) / total_surface_area,
    where *a_i* and *v_i* are the area and value for each vertex *i*.

    Parameters
    ----------
    values : numpy array of integer or float values
    areas : numpy array of surface areas
    labels : array or list of integer labels (same length as values)
    exclude_values : list of integer labels to be excluded

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
    >>> from mindboggle.utils.io_vtk import load_scalars
    >>> from mindboggle.measure.measure_functions import mean_value_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                            'labels', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> points, faces, areas, n_vertices = load_scalars(area_file, True)
    >>> points, faces, labels, n_vertices = load_scalars(label_file, True)
    >>> exclude_values = [-1,0]
    >>> mean_values, norm_mean_values, surface_areas, \
    >>>     label_list = mean_value_per_label(depths, areas, labels, exclude_values)

    """
    import numpy as np

    def avg_by_area(values_label, areas_label):
        return sum(areas_label * values_label) / sum(areas_label)

    label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_values]
    mean_values = []
    norm_mean_values = []
    surface_areas = []

    for label in label_list:

        I = [i for i,x in enumerate(labels) if x == label]
        mean_value = np.mean(values[I])
        norm_mean_value = avg_by_area(values[I], areas[I])
        mean_values.append(mean_value)
        norm_mean_values.append(norm_mean_value)
        surface_area = sum(areas[I])
        surface_areas.append(surface_area)

    return mean_values, norm_mean_values, surface_areas, label_list

def compute_percentile(N, percent, key=lambda x:x):
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
