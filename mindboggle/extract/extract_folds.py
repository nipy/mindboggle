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
#from measure.measure_functions import compute_percentile
#from utils.mesh_operations import segment, fill_holes
#from utils.io_vtk import load_scalar

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
    >>> from utils.mesh_operations import find_all_neighbors_from_file
    >>> from extract.extract_folds import extract_folds
    >>> depth_file = '/desk/output/results/measures/_hemi_rh_subject_MMRR-21-1/rh.pial.depth.vtk'
    >>> neighbor_lists = find_all_neighbors_from_file(depth_file)
    >>> folds, n_folds = extract_folds(depth_file, neighbor_lists, 0.5, 50)
    >>> # Write results to vtk file:
    >>> from utils.io_vtk import rewrite_scalars
    >>> rewrite_scalars(depth_file, 'test_folds.vtk', folds)

    """
    import numpy as np
    from time import time
    from measure.measure_functions import compute_percentile
    from utils.mesh_operations import segment, fill_holes
    from utils.io_vtk import load_scalar

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
    seeds = np.where(depths > min_depth)[0]
    folds, n_folds, max_fold, = segment(seeds, neighbor_lists, min_fold_size)
    print('    ...Folds segmented ({0:.2f} seconds)'.format(time() - t1))

    # If there are any folds
    if n_folds > 0:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep) and have some minimum depth
        t2 = time()
        seeds = [i for i,x in enumerate(folds) if x==0 if depths[i] > min_depth]

        # Segment holes in the folds
        print('  Segment holes in the folds...')
        holes, n_holes, max_hole = segment(seeds, neighbor_lists, 1)

        # If there are any holes
        if n_holes > 0:

            # Ignore the largest hole (the background) and renumber holes
            holes[holes == max_hole] = 0
            if max_hole < n_holes:
                holes[holes > max_hole] -= 1
            n_holes -= 1
            print('    ...Holes segmented ({0:.2f} seconds)'.format(time() - t2))

            t3 = time()
            folds = fill_holes(folds, holes, n_holes, neighbor_lists)
            print('  Filled holes ({0:.2f} seconds)'.format(time() - t3))

    print('  ...Extracted folds greater than {0:.2f} depth in {1:.2f} seconds'.
          format(min_depth, time() - t0))

    # Return folds, number of folds
    return folds, n_folds
