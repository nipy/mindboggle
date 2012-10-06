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

#==============
# Extract folds
#==============
def extract_folds(depth_file, neighbor_lists, fraction_folds, min_fold_size):
    """
    Use depth to extract folds from a triangular surface mesh and fill holes
    resulting from shallower areas within a fold.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    fraction_folds : float
        fraction of surface mesh considered folds
    min_fold_size : int
        minimum fold size (number of vertices)

    Returns
    -------
    folds : numpy array
        fold IDs
    n_folds :  int
        number of folds

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, write_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors_from_file
    >>> from mindboggle.extract.extract_folds import extract_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=0)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> folds, n_folds = extract_folds(depth_file, neighbor_lists, 0.5, 50)
    >>> # Write results to vtk file:
    >>> from mindboggle.utils.io_vtk import rewrite_scalars
    >>> rewrite_scalars(depth_file, 'test_extract_folds.vtk', folds)

    """
    import numpy as np
    from time import time
    from mindboggle.measure.measure_functions import compute_percentile
    from mindboggle.utils.mesh_operations import segment, fill_holes
    from mindboggle.utils.io_vtk import load_scalar

    print("Extract folds from surface mesh...")
    t0 = time()

    # Load depth values from VTK file
    points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=True)

    # Compute the minimum depth threshold for defining folds by determining the
    # percentile of depth values for the fraction of vertices that are not folds.
    # For example, if we consider the shallowest one-third of vertices not to be
    # folds, we compute the depth percentile, and two-thirds of vertices would
    # have at least this depth value and would be considered folds.
    min_depth = compute_percentile(np.sort(depths),
                                   1 - fraction_folds, key=lambda x:x)

    # Segment folds of a surface mesh
    print("  Segment surface mesh into separate folds deeper than {0:.2f}...".
          format(min_depth))
    t1 = time()
    vertices_to_segment = np.where(depths > min_depth)[0]
    folds, n_folds = segment(vertices_to_segment, neighbor_lists,
        seed_lists=[], min_region_size=min_fold_size,
        spread_same_labels=False, labels=[], label_pair_lists=[])
    print('    ...Folds segmented ({0:.2f} seconds)'.format(time() - t1))

    # If there are any folds
    if n_folds > 0:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep)
        t2 = time()
        vertices_to_segment = [i for i,x in enumerate(folds) if x==-1]

        # Segment holes in the folds
        print('  Segment holes in the folds...')
        holes, n_holes = segment(vertices_to_segment, neighbor_lists,
                seed_lists=[], min_region_size=1,
                spread_same_labels=False, labels=[], label_pair_lists=[])

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
            folds = fill_holes(folds, holes, n_holes, neighbor_lists)
            print('  Filled holes ({0:.2f} seconds)'.format(time() - t3))

    print('  ...Extracted folds greater than {0:.2f} depth in {1:.2f} seconds'.
          format(min_depth, time() - t0))

    # Return folds, number of folds
    return folds, n_folds
