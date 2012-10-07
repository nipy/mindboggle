#!/usr/bin/python
"""
Function to extract folds from a triangular surface mesh.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import numpy as np
#from time import time
#from mindboggle.measure.measure_functions import compute_percentile
#from mindboggle.utils.mesh_operations import segment, fill_holes
#from mindboggle.utils.io_vtk import load_scalar

#===============================================================================
# Extract folds
#===============================================================================
def extract_all_folds(depth_file, area_file, fraction_folds):
    """
    Use depth to extract all folds from a triangular surface mesh,
    by fraction of surface area.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    area_file : str
        surface mesh file in VTK format with faces and surface area scalar values
    fraction_folds : float
        fraction of surface mesh considered folds

    Returns
    -------
    folds : array of integers
        an integer for every mesh vertex: 1 for fold, -1 for non-fold

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, write_scalars
    >>> from mindboggle.extract.extract_folds import extract_all_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> folds, min_depth, n_deep_vertices = extract_all_folds(depth_file,
    >>>     area_file, 0.5)
    >>> # Write results to vtk file and view with mayavi2:
    >>> from mindboggle.utils.io_vtk import rewrite_scalars
    >>> rewrite_scalars(depth_file, 'test_extract_all_folds.vtk', folds, folds)
    >>> os.system('mayavi2 -m Surface -d test_extract_all_folds.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.utils.io_vtk import load_scalar

    print("Extract the deepest surface mesh vertices ({0} of surface area)...".
          format(fraction_folds))
    t0 = time()

    # Load depth and surface area values from VTK files
    points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=True)
    points, faces, areas, n_vertices = load_scalar(area_file, return_arrays=True)

    sort_indices = np.argsort(depths)
    sort_indices = sort_indices[::-1]

    total_area = np.sum(areas)
    fraction_area = fraction_folds * total_area
    sum_area = 0
    folds = -1 * np.ones(len(areas))
    min_depth = 1
    for index in sort_indices:
        folds[index] = 1
        sum_area += areas[index]
        if sum_area >= fraction_area:
            min_depth = depths[index]
            break

    n_deep_vertices = len([x for x in folds if x == 1])

    print('  ...Extracted {0} vertices deeper than {1:.2f} ({2:.2f} seconds)'.
          format(n_deep_vertices, min_depth, time() - t0))

    return folds, min_depth, n_deep_vertices

#===============================================================================
# Extract individual folds
#===============================================================================
def extract_folds(depth_file, area_file, neighbor_lists, fraction_folds, min_fold_size):
    """
    Use depth to extract folds from a triangular surface mesh and fill holes
    resulting from shallower areas within a fold.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    area_file : str
        surface mesh file in VTK format with faces and surface area scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    fraction_folds : float
        fraction of surface mesh considered folds
    min_fold_size : int
        minimum fold size (number of vertices)

    Returns
    -------
    fold_IDs : array of integers
        fold IDs for all vertices, with -1s for non-fundus vertices
    n_folds :  int
        number of folds

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.extract.extract_folds import extract_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=0)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> fold_IDs, n_folds = extract_folds(depth_file, area_file, neighbor_lists, 0.5, 50)
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalars(depth_file, 'test_extract_folds.vtk', fold_IDs, fold_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_folds.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.measure.measure_functions import compute_percentile
    from mindboggle.utils.mesh_operations import segment, fill_holes
    from mindboggle.utils.io_vtk import load_scalar
    from mindboggle.extract.extract_folds import extract_all_folds

    print("Extract folds from surface mesh...")
    t0 = time()

    # Compute the minimum depth threshold for defining folds
    folds, min_depth, n_deep_vertices = extract_all_folds(depth_file,
        area_file, fraction_folds)

    # Segment folds of a surface mesh
    print("  Segment surface mesh into separate folds deeper than {0:.2f}...".
          format(min_depth))
    t1 = time()
    vertices_to_segment = [i for i,x in enumerate(folds) if x == 1]
    fold_IDs = segment(vertices_to_segment, neighbor_lists,
        seed_lists=[], min_region_size=min_fold_size,
        spread_same_labels=False, labels=[], label_pair_lists=[])
    print('    ...Folds segmented ({0:.2f} seconds)'.format(time() - t1))
    n_folds = len([x for x in list(set(fold_IDs))])

    # If there are any folds
    if n_folds > 0:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep)
        t2 = time()
        vertices_to_segment = [i for i,x in enumerate(fold_IDs) if x==-1]

        # Segment holes in the folds
        print('  Segment holes in the folds...')
        holes = segment(vertices_to_segment, neighbor_lists,
                seed_lists=[], min_region_size=1,
                spread_same_labels=False, labels=[], label_pair_lists=[])
        n_holes = len([x for x in list(set(holes))])

        # If there are any holes
        if n_holes > 0:

            # Ignore the largest hole (the background) and renumber holes
            max_hole_size = 0
            max_hole_index = 0
            for ihole in range(n_holes):
                I = np.where(holes == ihole)
                if len(I) > max_hole_size:
                    max_hole_size = len(I)
                    max_hole_index = ihole
            holes[holes == max_hole_index] = -1
            if max_hole_index < n_holes:
                holes[holes > max_hole_index] -= 1
            n_holes -= 1
            print('    ...{0} holes segmented ({1:.2f} seconds)'.
                  format(n_holes, time() - t2))

            # Fill holes
            t3 = time()
            fold_IDs = fill_holes(fold_IDs, holes, n_holes, neighbor_lists)
            print('  Filled holes ({0:.2f} seconds)'.format(time() - t3))

    print('  ...Extracted folds greater than {0:.2f} depth in {1:.2f} seconds'.
          format(min_depth, time() - t0))

    # Return folds, number of folds
    return fold_IDs, n_folds
