#!/usr/bin/env python
"""
Shape calculations.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

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
    Example: CurvatureMain -m 2 -n 0.7  lh.pial.vtk  lh.pial.mean_curvature.vtk
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
    >>> from mindboggle.shapes.measure import means_per_label
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
                label_weight = sum(W)
                label_areas.append(label_weight)
                if dim > 1:
                    W = np.transpose(np.tile(W, (dim,1)))
                means.append(np.sum(W * X, axis=0) / label_weight)
                sdevs.append(np.sqrt(np.sum(W * (X - np.mean(X, axis=0))**2,
                                            axis=0) / label_weight))
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
    >>> from mindboggle.shapes.measure import stats_per_label
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
    >>> #medians, mads, means, sdevs, skews, kurts, \
    >>> #lower_quarts, upper_quarts, label_list  = stats_per_label(values,
    >>> #    labels, exclude_labels, weights, precision)
    >>> stats_per_label(values, labels, exclude_labels, weights, precision)

    """
    import numpy as np
    from scipy.stats import skew, kurtosis, scoreatpercentile
    from mindboggle.utils.compute import weighted_to_repeated_values, median_abs_dev

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)
    if not isinstance(weights, np.ndarray):
        weights = np.asarray(weights)

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

    for label in label_list:
        I = [i for i,x in enumerate(labels) if x == label]
        if I:
            X = values[I]
            if np.size(weights):
                W = weights[I]
                sumW = np.sum(W)
                Xdiff = X - np.mean(X)
                means.append(np.sum(W * X) / sumW)
                sdevs.append(np.sum(W * Xdiff**2) / sumW)
                skews.append(np.sum(W * Xdiff**3) / sumW)
                kurts.append(np.sum(W * Xdiff**4) / sumW)
                X = weighted_to_repeated_values(X, weights[I], precision)
            else:
                means.append(np.mean(X))
                sdevs.append(np.std(X))
                skews.append(skew(X))
                kurts.append(kurtosis(X))
            medians.append(np.median(X))
            mads.append(median_abs_dev(X))
            lower_quarts.append(scoreatpercentile(X, 25))
            upper_quarts.append(scoreatpercentile(X, 75))
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
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
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
    >>> input_vtk = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
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
