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

def travel_depth(command, surface_file):
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
                 os.path.splitext(os.path.basename(surface_file))[0] + '.travel_depth.vtk')
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, depth_file])
    cli.cmdline
    cli.run()

    return depth_file

def geodesic_depth(command, surface_file):
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
                 os.path.splitext(os.path.basename(surface_file))[0] + '.geodesic_depth.vtk')
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, depth_file])
    cli.cmdline
    cli.run()

    return depth_file

def curvature(command, method, arguments, surface_file):
    """
    Measure curvature values of each vertex in a surface mesh (-m 0).
    (Calls Joachim Giard's C++ code)

    Usage: CurvatureMain [Options] InputVTKMesh MeanCurvatureOutput
    Options:
       -m Method: set the method used to compute curvature(s) (default 0)
           0 -- Use ComputePrincipalCurvatures() function to compute both
                mean and Gaussian curvatures based on the relative direction
                of the normal vectors in a small neighborhood
           1 -- Use ComputeBothCurvatures() function to compute both mean
                and Gaussian curvatures based on the local ratios between
                a filtered surface and the original surface area
           2 -- Use ComputeCurvature() function to compute the mean curvature
                based on the direction of the displacement vectors during
                a Laplacian filtering
       -n Neighborhood: neighborhood size (default 0.7)
       -g GaussianCurvVTK: save Gaussian curvature (for -m 0 or 1)
       -x MaxCurvVTK: save maximum curvature (for -m 0)
       -i MinCurvVTK: save minimum curvature (for -m 0)
       -d DirectionVTK: save minimal curvature's direction (for -m 0)
    Example: CurvatureMain -m 2 -n 0.7  lh.pial.vtk  lh.mean_curv.vtk
    Example: CurvatureMain -m 0 -n 2
                -i lh.min_curv.vtk -x lh.max_curv.vtk -g lh.gaussian_curv.vtk
                -d lh.min_dir.vtk lh.pial.vtk  lh.mean_curv.vtk

    Note ::

        -m 0 is best if you have low resolution or want local peaks,
            but can be too sensitive to the local linear geometry of the mesh,
            unless the neighborhood parameter is set high enough (like 2).
        -m 1 is not well tested and the filtering is done using Euclidean
            distances, so it's only good for incorrect but fast visualization.
        -m 2 is a good approximation based on the Laplacian, but very large
            curvatures (negative or positive) are underestimated (saturation).

    Parameters
    ----------
    command : string
        C++ executable command for computing curvature
    method : integer {0,1,2}
        method number
    arguments : string
        additional arguments, such as neighborhood parameter
    surface_file : string
        name of VTK surface mesh file

    """
    import os
    from nipype.interfaces.base import CommandLine

    args = ['-m', str(method)]
    gauss_curvature_file = None
    max_curvature_file = None
    min_curvature_file = None
    min_curvature_vector_file = None

    stem = os.path.join(os.getcwd(),
                        os.path.splitext(os.path.basename(surface_file))[0])
    mean_curvature_file = stem + '.mean_curvature.vtk'
    if method in [0,1]:
        gauss_curvature_file = stem + '.gauss_curvature.vtk'
        args.extend(['-g', gauss_curvature_file])
    if method == 0:
        max_curvature_file = stem + '.max_curvature.vtk'
        min_curvature_file = stem + '.min_curvature.vtk'
        min_curvature_vector_file = stem + '.min_curvature.txt'
        args.extend(['-x', max_curvature_file,
                     '-i', min_curvature_file,
                     '-d', min_curvature_vector_file])
    args.extend([surface_file, mean_curvature_file])

    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    return mean_curvature_file, gauss_curvature_file, \
           max_curvature_file, min_curvature_file, min_curvature_vector_file

def mean_value_per_label(values, labels, exclude_labels,
                         normalize_by_area=False, areas=[]):
    """
    Compute the mean value across vertices per label,
    optionally taking into account surface area per vertex.

    average value = sum(a_i * v_i) / total_surface_area,
    where *a_i* and *v_i* are the area and value for each vertex *i*.

    Parameters
    ----------
    values : numpy array of individual or lists of integers or floats
        values to average per label
    labels : list or array of integers
        label for each value
    exclude_labels : list of integers
        labels to be excluded
    normalize_by_area : Boolean
        divide each mean value per label by the surface area of that label?
    areas : numpy array of floats (if normalize_by_area)
        surface areas

    Returns
    -------
    mean_values : list of floats
        mean value(s) for each label
    label_list : list of integers
        unique label numbers
    label_areas : list of floats (if normalize_by_area)
        surface area for each labeled set of vertices
    norm_mean_values : list of floats (if normalize_by_area)
        mean values normalized by vertex area

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.shapes.measure import mean_value_per_label
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> labels_file = os.path.join(data_path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> values, name = read_scalars(depth_file, True, True)
    >>> areas, name = read_scalars(area_file, True, True)
    >>> labels, name = read_scalars(labels_file)
    >>> exclude_labels = [-1]
    >>> normalize_by_area = True
    >>> mean_values, label_list, label_areas, norm_mean_values = mean_value_per_label(values,
    >>>     labels, exclude_labels, normalize_by_area, areas)

    """
    import numpy as np

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(areas, np.ndarray):
        areas = np.asarray(areas)

    label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) not in exclude_labels]
    mean_values = []
    label_areas = []
    norm_mean_values = []
    if values.ndim > 1:
        dim = np.shape(values)[1]
    else:
        dim = 1

    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            if dim > 1:
                mean_value = np.mean(values[I], axis=0)
            else:
                mean_value = np.mean(values[I])
            mean_values.append(mean_value)

            if normalize_by_area:
                surface_area = sum(areas[I])
                label_areas.append(surface_area)
                if dim > 1:
                    areas_dup = np.transpose(np.tile(areas[I], (dim,1)))
                    norm_mean_value = sum(values[I] * areas_dup) / sum(areas[I])
                else:
                    norm_mean_value = sum(values[I] * areas[I]) / sum(areas[I])
                norm_mean_values.append(norm_mean_value)
        else:
            mean_values.append(0)
            if normalize_by_area:
                label_areas.append(0)
                norm_mean_values.append(0)

    mean_values = [x.tolist() for x in mean_values]
    label_areas = [x.tolist() for x in label_areas]
    norm_mean_values = [x.tolist() for x in norm_mean_values]

    return mean_values, label_list, label_areas, norm_mean_values

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
    >>> labels = [int(x) for x in labels]
    >>> volumes, labels = volume_per_label(labels, input_file)
    >>> print(volumes)

    """
    import numpy as np
    import nibabel as nb

    # Load labeled image volumes
    img = nb.load(input_file)
    hdr = img.get_header()
    pixdims = hdr.get_zooms()
    volume_per_voxel = np.product(pixdims)
    data = img.get_data().ravel()

    # Initialize output
    volumes = np.zeros(len(labels))

    # Loop through labels
    for ilabel, label in enumerate(labels):
        label = int(label)

        # Find which voxels contain the label in each volume
        indices = np.where(data==label)[0]
        volumes[ilabel] = volume_per_voxel * len(indices)

    return volumes.tolist(), labels

def rescale_by_neighborhood(input_vtk, indices=[], nedges=10, p=99,
    set_max_to_1=True, save_file=False, output_filestring='rescaled_scalars'):
    """
    Rescale the scalar values of a VTK file by a percentile value
    in each vertex's surface mesh neighborhood.

    Parameters
    ----------
    input_vtk : string
        name of VTK file with a scalar value for each vertex
    indices : list of integers (optional)
        indices of scalars to normalize
    nedges : integer
        number or edges from vertex, defining the size of its neighborhood
    p : float in range of [0,100]
        percentile used to normalize each scalar
    set_max_to_1 : Boolean
        set all rescaled values greater than 1 to 1.0?
    save_file : Boolean
        save output VTK file?
    output_filestring : string (if save_file)
        name of output file

    Returns
    -------
    rescaled_scalars : list of floats
        rescaled scalar values
    rescaled_scalars_file : string (if save_file)
        name of output VTK file with rescaled scalar values

    Examples
    --------
    >>> import os
    >>> from mindboggle.shapes.measure import rescale_by_neighborhood
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> indices = []
    >>> nedges = 10
    >>> p = 99
    >>> set_max_to_1 = True
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> #
    >>> rescaled_scalars, rescaled_scalars_file = rescale_by_neighborhood(input_vtk,
    >>>     indices, nedges, p, set_max_to_1, save_file, output_filestring)
    >>> #
    >>> # View rescaled scalar values per fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> rewrite_scalars(rescaled_scalars_file, rescaled_scalars_file,
    >>>                 rescaled_scalars, 'rescaled_depths', folds)
    >>> plot_vtk(rescaled_scalars_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file, find_neighborhood
    from mindboggle.shapes.measure import rescale_by_neighborhood

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    if not indices:
        indices = [i for i,x in enumerate(scalars) if x != -1]
    print("  Rescaling {0} scalar values by neighborhood...".format(len(indices)))
    neighbor_lists = find_neighbors_from_file(input_vtk)

    # Loop through vertices:
    rescaled_scalars = scalars.copy()
    for index in indices:

        # Determine the scalars in the vertex's neighborhood:
        neighborhood = find_neighborhood(neighbor_lists, [index], nedges)

        # Compute a high neighborhood percentile to normalize vertex's value:
        normalization_factor = np.percentile(scalars[neighborhood], p)
        rescaled_scalar = scalars[index] / normalization_factor
        rescaled_scalars[index] = rescaled_scalar

    # Make any rescaled value greater than 1 equal to 1:
    if set_max_to_1:
        rescaled_scalars[[x for x in indices if rescaled_scalars[x] > 1.0]] = 1

    rescaled_scalars = rescaled_scalars.tolist()

    #---------------------------------------------------------------------------
    # Return rescaled scalars and file name
    #---------------------------------------------------------------------------
    if save_file:

        rescaled_scalars_file = os.path.join(os.getcwd(), output_filestring + '.vtk')
        rewrite_scalars(input_vtk, rescaled_scalars_file,
                        rescaled_scalars, 'rescaled_scalars')

    else:
        rescaled_scalars_file = None

    return rescaled_scalars, rescaled_scalars_file


def rescale_by_label(input_vtk, labels_or_file, combine_all_labels=False,
                     nedges=10, p=99, set_max_to_1=True, save_file=False,
                     output_filestring='rescaled_scalars'):
    """
    Rescale scalars for each label (such as depth values within each fold).

    Default is to normalize the scalar values of a VTK file by
    a percentile value in each vertex's surface mesh for each label.

    Parameters
    ----------
    input_vtk : string
        name of VTK file with a scalar value for each vertex
    labels_or_file : list or string
        label number for each vertex or name of VTK file with index scalars
    combine_all_labels : Boolean
        combine all labels (scalars not equal to -1) as one label?
    nedges : integer (if norm_by_neighborhood)
        number or edges from vertex, defining the size of its neighborhood
    p : float in range of [0,100] (if norm_by_neighborhood)
        percentile used to rescale each scalar
    set_max_to_1 : Boolean
        set all rescaled values greater than 1 to 1.0?
    save_file : Boolean
        save output VTK file?
    output_filestring : string (if save_file)
        name of output file

    Returns
    -------
    rescaled_scalars : list of floats
        scalar values rescaled for each label, for label numbers not equal to -1
    rescaled_scalars_file : string (if save_file)
        name of output VTK file with rescaled scalar values for each label

    Examples
    --------
    >>> # Rescale depths by neighborhood within each label:
    >>> import os
    >>> from mindboggle.shapes.measure import rescale_by_label
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> labels_or_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> combine_all_labels = False
    >>> nedges = 10
    >>> p = 99
    >>> set_max_to_1 = True
    >>> save_file = True
    >>> output_filestring = 'rescaled_scalars'
    >>> #
    >>> rescaled_scalars, rescaled_scalars_file = rescale_by_label(input_vtk,
    >>>     labels_or_file, combine_all_labels, nedges, p,
    >>>     set_max_to_1, save_file, output_filestring)
    >>> #
    >>> # View rescaled scalar values per fold:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> rewrite_scalars(rescaled_scalars_file, rescaled_scalars_file,
    >>>                 rescaled_scalars, 'rescaled_depths', folds)
    >>> plot_vtk(rescaled_scalars_file)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file

    # Load scalars and vertex neighbor lists:
    scalars, name = read_scalars(input_vtk, True, True)
    print("  Rescaling scalar values within each label...")

    # Load label numbers:
    if isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file, True, True)
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    unique_labels = np.unique(labels)
    unique_labels = [x for x in unique_labels if x >= 0]

    # Loop through labels:
    for label in unique_labels:
        #print("  Rescaling scalar values within label {0} of {1} labels...".format(
        #    int(label), len(unique_labels)))
        indices = [i for i,x in enumerate(labels) if x == label]
        if indices:

            # Rescale by the maximum label scalar value:
            scalars[indices] = scalars[indices] / np.max(scalars[indices])

    rescaled_scalars = scalars.tolist()

    #---------------------------------------------------------------------------
    # Return rescaled scalars and file name
    #---------------------------------------------------------------------------
    if save_file:

        rescaled_scalars_file = os.path.join(os.getcwd(), output_filestring + '.vtk')
        rewrite_scalars(input_vtk, rescaled_scalars_file,
                        rescaled_scalars, 'rescaled_scalars', labels)

    else:
        rescaled_scalars_file = None

    return rescaled_scalars, rescaled_scalars_file

