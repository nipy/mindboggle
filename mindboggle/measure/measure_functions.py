#!/usr/bin/python

"""
Shape calculations.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

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

def mean_value_per_label(values, areas, labels):
    """
    Compute the mean value across vertices per label,
    taking into account surface area per vertex.

    Parameters
    ----------
    values : numpy array of integer or float values
    areas : numpy array of surface areas
    labels : array or list of integer labels (same length as values)

    Returns
    -------
    mean_values : list of floats
    label_list : list of unique labels

    """
    import numpy as np

    # Weight each value by the surface area for the given vertex
    total_surface_area = sum(areas)
    values = values * areas / (len(values) * total_surface_area)

    label_list = np.unique(labels)
    label_list = [int(x) for x in label_list if int(x) != 0]
    mean_values = []
    surface_areas = []

    for label in label_list:

        I = [i for i,x in enumerate(labels) if x == label]
        mean_value = np.mean(values[I])
        mean_values.append(mean_value)
        surface_area = sum(areas[I])
        surface_areas.append(surface_area)

    return mean_values, surface_areas, label_list
