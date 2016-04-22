#!/usr/bin/env python
"""
Shape calculations.


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def area(command, surface_file, verbose=False):
    """
    Measure area of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    Parameters
    ----------
    command : string
        Voronoi-based surface area C++ executable command
    surface_file : string
        vtk file with surface mesh
    verbose : bool
        print statements?

    Returns
    -------
    area_file: string
        vtk file with surface area per vertex of mesh

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.surface_shapes import area
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> surface_file = fetch_data(urls['left_pial'])
    >>> verbose = False
    >>> ccode_path = os.environ['vtk_cpp_tools']
    >>> command = os.path.join(ccode_path, 'area', 'PointAreaMain')
    >>> area_file = area(command, surface_file, verbose)
    >>> scalars, name = read_scalars(area_file)
    >>> print(np.array_str(np.array(scalars[0:8]), precision=5,
    ...     suppress_small=True))
    [ 0.4827   0.39661  0.57813  0.70574  0.84318  0.57643  0.66942  0.7063 ]

    """
    import os
    from nipype.interfaces.base import CommandLine

    basename = os.path.splitext(os.path.basename(surface_file))[0]
    area_file = os.path.join(os.getcwd(), basename + '.area.vtk')
    args = ' '.join([surface_file, area_file])

    if verbose:
        print("{0} {1}".format(command, args))

    cli = CommandLine(command=command)
    cli.inputs.args = args
    cli.cmdline
    cli.run()

    if not os.path.exists(area_file):
        raise IOError(area_file + " not found")

    return area_file


def travel_depth(command, surface_file, verbose=False):
    """
    Measure "travel depth" of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    Parameters
    ----------
    command : string
        travel depth C++ executable command
    surface_file : string
        vtk file
    verbose : bool
        print statements?

    Returns
    -------
    depth_file: string
        vtk file with travel depth per vertex of mesh

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.surface_shapes import travel_depth
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> surface_file = fetch_data(urls['left_pial'])
    >>> verbose = False
    >>> ccode_path = os.environ['vtk_cpp_tools']
    >>> command = os.path.join(ccode_path, 'travel_depth', 'TravelDepthMain')
    >>> depth_file = travel_depth(command, surface_file, verbose)
    >>> scalars, name = read_scalars(depth_file)
    >>> print(np.array_str(np.array(scalars[0:8]), precision=5,
    ...     suppress_small=True))
    [ 0.02026  0.06009  0.12859  0.04564  0.00774  0.05284  0.05354  0.01316]

    """
    import os
    from nipype.interfaces.base import CommandLine

    basename = os.path.splitext(os.path.basename(surface_file))[0]
    depth_file = os.path.join(os.getcwd(), basename + '.travel_depth.vtk')
    args = ' '.join([surface_file, depth_file])

    if verbose:
        print("{0} {1}".format(command, args))

    cli = CommandLine(command=command)
    cli.inputs.args = args
    cli.cmdline
    cli.run()

    if not os.path.exists(depth_file):
        raise IOError(depth_file + " not found")

    return depth_file


def geodesic_depth(command, surface_file, verbose=False):
    """
    Measure "travel depth" of each vertex in a surface mesh.
    (Calls Joachim Giard's C++ code)

    Parameters
    ----------
    command : travel depth C++ executable command
    surface_file : ``vtk file``
    verbose : bool
        print statements?

    Returns
    -------
    depth_file: string
        vtk file with geodesic depth per vertex of mesh

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.surface_shapes import geodesic_depth
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> surface_file = fetch_data(urls['left_pial'])
    >>> verbose = False
    >>> ccode_path = os.environ['vtk_cpp_tools']
    >>> command = os.path.join(ccode_path, 'geodesic_depth', 'GeodesicDepthMain')
    >>> depth_file = geodesic_depth(command, surface_file, verbose)
    >>> scalars, name = read_scalars(depth_file)
    >>> print(np.array_str(np.array(scalars[0:8]), precision=5,
    ...     suppress_small=True))
    [ 0.02026  0.06009  0.12859  0.04564  0.00774  0.05284  0.05354  0.01316]

    """
    import os
    from nipype.interfaces.base import CommandLine

    basename = os.path.splitext(os.path.basename(surface_file))[0]
    depth_file = os.path.join(os.getcwd(), basename + '.geodesic_depth.vtk')
    args = ' '.join([surface_file, depth_file])

    if verbose:
        print("{0} {1}".format(command, args))

    cli = CommandLine(command=command)
    cli.inputs.args = args
    cli.cmdline
    cli.run()

    if not os.path.exists(depth_file):
        raise IOError(depth_file + " not found")

    return depth_file


def curvature(command, method, arguments, surface_file, verbose=False):
    """
    Measure curvature values of each vertex in a surface mesh (-m 0).
    (Calls Joachim Giard's C++ code)

    Command line usage:
    CurvatureMain [Options] InputVTKMesh MeanCurvatureOutput

    Command line options:
    -m Method: set the method used to compute curvature(s) (default 0):
    0: Use ComputePrincipalCurvatures() function to compute both
    mean and Gaussian curvatures based on the relative direction
    of the normal vectors in a small neighborhood
    1: Use ComputeBothCurvatures() function to compute both mean
    and Gaussian curvatures based on the local ratios between
    a filtered surface and the original surface area
    2: Use ComputeCurvature() function to compute the mean curvature
    based on the direction of the displacement vectors during
    a Laplacian filtering
    -n Neighborhood: neighborhood size (default 0.7)
    -g GaussianCurvVTK: save Gaussian curvature (for -m 0 or 1)
    -x MaxCurvVTK: save maximum curvature (for -m 0)
    -i MinCurvVTK: save minimum curvature (for -m 0)
    -d DirectionVTK: save minimal curvature's direction (for -m 0)

    Command line example:
    CurvatureMain -m 2 -n 0.7  lh.pial.vtk  lh.pial.mean_curvature.vtk

    Command line example:
    CurvatureMain -m 0 -n 2
    -i lh.min_curv.vtk -x lh.max_curv.vtk -g lh.gaussian_curv.vtk
    -d lh.min_dir.vtk lh.pial.vtk  lh.mean_curv.vtk

    Notes:
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
    verbose : bool
        print statements?

    Returns
    -------
    mean_curvature_file : string
        mean curvature file
    gauss_curvature_file : string
        gauss curvature file
    max_curvature_file : string
        max curvature file
    min_curvature_file : string
        min curvature file
    min_curvature_vector_file : string
        min curvature vector file

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.surface_shapes import curvature
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> surface_file = fetch_data(urls['left_pial'])
    >>> method = 2
    >>> arguments = '-n 0.7'
    >>> verbose = False
    >>> ccode_path = os.environ['vtk_cpp_tools']
    >>> command = os.path.join(ccode_path, 'curvature', 'CurvatureMain')
    >>> mean_curvature_file, f1,f2,f3,f4 = curvature(command, method,
    ...     arguments, surface_file, verbose)
    >>> scalars, name = read_scalars(mean_curvature_file)
    >>> print(np.array_str(np.array(scalars[0:8]), precision=5,
    ...     suppress_small=True))
    [-5.81361 -5.9313  -6.28055 -5.621   -5.69631 -5.80399 -5.87265 -5.7107 ]

    """
    import os
    from nipype.interfaces.base import CommandLine

    args = ['-m', str(method)]
    gauss_curvature_file = None
    max_curvature_file = None
    min_curvature_file = None
    min_curvature_vector_file = None

    basename = os.path.splitext(os.path.basename(surface_file))[0]
    mean_curvature_file = os.path.join(os.getcwd(), basename) + \
        '.mean_curvature.vtk'
    if method in [0, 1]:
        gauss_curvature_file = stem + '.gauss_curvature.vtk'
        args.extend(['-g', gauss_curvature_file])
    if method == 0:
        max_curvature_file = stem + '.max_curvature.vtk'
        min_curvature_file = stem + '.min_curvature.vtk'
        min_curvature_vector_file = stem + '.min_curvature.txt'
        args.extend(['-x', max_curvature_file,
                     '-i', min_curvature_file,
                     '-d', min_curvature_vector_file])

    if arguments:
        args.extend([arguments])

    args.extend([surface_file, mean_curvature_file])

    if verbose:
        print("{0} {1}".format(command, args))

    cli = CommandLine(command=command)
    cli.inputs.args = ' '.join(args)
    cli.cmdline
    cli.run()

    return mean_curvature_file, gauss_curvature_file, \
           max_curvature_file, min_curvature_file, min_curvature_vector_file


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules