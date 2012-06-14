#!/usr/bin/python

"""
Surface calculations


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

##############################################################################
#   Surface calculations
##############################################################################

def compute_depth(command, surface_file):
    """
    Measure

    measure_()
    """
    from os import getcwd, path, error
    from nipype.interfaces.base import CommandLine

    # Check type:
    if type(surface_file) == str:
        pass
    elif type(surface_file) == list:
        surface_file = surface_file[0]
    else:
        error("Check format of " + surface_file)

    depth_file = path.splitext(path.basename(surface_file))[0] + '.depth.vtk'
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([surface_file, path.join(getcwd(), depth_file)])
    cli.cmdline
    cli.run()
    return depth_file
    
def compute_curvature(command, surface_file):
    """
    Measure
    CurvatureMain input MeanCurvatureOutput [GaussianCurvatureOutput]
                        [MaximalCurvatureOutput] [MinimalCurvatureOutput]
    measure_()
    """
    from os import getcwd, path, error
    from nipype.interfaces.base import CommandLine

    # Check type:
    if type(surface_file) == str:
        pass
    elif type(surface_file) == list:
        surface_file = surface_file[0]
    else:
        error("Check format of " + surface_file)

    file_stem = path.join(getcwd(), path.splitext(path.basename(surface_file))[0])
    mean_curvature_file = file_stem + '.curvature.mean.vtk'
    gauss_curvature_file = file_stem + '.curvature.gauss.vtk'
    max_curvature_file = file_stem + '.curvature.max.vtk'
    min_curvature_file = file_stem + '.curvature.min.vtk'
    min_curvature_vector_file = file_stem + '.curvature.min.direction.txt'
    args = [surface_file,
            mean_curvature_file, gauss_curvature_file,
            max_curvature_file, min_curvature_file]
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()
    return mean_curvature_file, gauss_curvature_file,\
           max_curvature_file, min_curvature_file, min_curvature_vector_file
