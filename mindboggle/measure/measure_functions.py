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
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
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
