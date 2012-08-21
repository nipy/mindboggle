#!/usr/bin/python

"""
Shape calculations


Authors:
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com
Forrest Sheng Bao  .  http://fsbao.net

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Surface calculations
##############################################################################

def compute_depth(command, surface_file):
    """
    Measure Joachim Giard's "travel depth" for a surface mesh.
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

def average_value_per_label(values, labels):
    """
    Compute the mean value per label.

    Inputs:
    ======
    values:  list of values
    labels:  list of integer labels (same length as values)

    Output:
    ======
    mean_values:  list of floats

    """
    import numpy as np

    label_list = np.unique(labels)
    mean_values = []
    for label in label_list:
        mean_value = np.mean(values[np.where(labels == label)[0]])
        mean_values.append(mean_value)

    return mean_values
