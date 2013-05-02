#!/usr/bin/env python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Arno Klein, 2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Extract fundi
#=============================================================================
def extract_fundi(folds_or_file, depth_file, likelihoods_or_file,
                  smooth_skeleton=False, save_file=False):
    """
    Extract fundi from folds.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Steps ::

        1. Find fundus endpoints using find_outer_anchors().
        2. Connect fundus endpoints using connect_points_erosion().

    Parameters
    ----------
    folds_or_file : list or string
        fold number for each vertex or name of VTK file containing folds scalars
    depth_file :  string
        surface mesh file in VTK format with scalar rescaled depth values
    likelihoods_or_file : list or string
        fundus likelihood values or name of VTK file with the scalar values
    smooth_skeleton : Boolean
        smooth skeleton?
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundi : list of integers
        fundus numbers for all vertices (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundus_endpoints : list of integers
        indices to fundus anchor vertices
    fundi_file : string (if save_file)
        name of output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.plots import plot_vtk
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> single_fold = True
    >>> if single_fold:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, True, True)
    >>>     fold_number = 11 #11
    >>>     folds[folds != fold_number] = -1
    >>> else:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> #
    >>> smooth_skeleton = False
    >>> fundi, n_fundi, endpoints, fundi_file = extract_fundi(folds, depth_file,
    >>>     likelihoods_or_file, smooth_skeleton, save_file=True)
    >>> #
    >>> # View:
    >>> plot_vtk(fundi_file)

    """
    import os
    import sys
    import numpy as np
    from time import time

    from mindboggle.utils.paths import find_outer_anchors, \
        connect_points_erosion, connect_points_hmmf
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.morph import dilate
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Load fold numbers:
    if isinstance(folds_or_file, str):
        folds, name = read_scalars(folds_or_file)
    elif isinstance(folds_or_file, list):
        folds = folds_or_file
    elif isinstance(folds_or_file, np.ndarray):
        folds = folds_or_file.tolist()
    else:
        sys.error('folds_or_file is not a string, list, or array.')

    # Load likelihood values:
    if isinstance(likelihoods_or_file, str):
        likelihoods, name = read_scalars(likelihoods_or_file, True, True)
    elif isinstance(likelihoods_or_file, list):
        likelihoods = np.array(likelihoods_or_file)
    elif isinstance(likelihoods_or_file, np.ndarray):
        likelihoods = likelihoods_or_file

    # Initialize variables:
    t1 = time()
    count = 0
    neighbor_lists = find_neighbors_from_file(depth_file)
    depths, name = read_scalars(depth_file)
    npoints = len(folds)
    Z = np.zeros(npoints)
    fundi = -1 * np.ones(npoints)
    fundus_endpoints = []
    nedges = 2

    # For each fold region...
    unique_fold_IDs = [x for x in np.unique(folds) if x != -1]
    print("Extract a fundus from each of {0} regions...".
          format(len(unique_fold_IDs)))
    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Fold {0}:'.format(int(fold_ID)))

            # Find fundus endpoints on the boundary of the surface region:
            endpoints, endtracks = find_outer_anchors(indices_fold,
                neighbor_lists, likelihoods, depths, min_edges=5)

            #-----------------------------------------------------------------
            # Connect endpoints to each other:
            #-----------------------------------------------------------------
            # Connect endpoints via a skeleton:
            B = -1 * np.ones(len(folds))
            B[indices_fold] = 1
            skeleton = connect_points_erosion(B, endpoints, neighbor_lists,
                                              likelihoods, test_ratio=0.5)
            #skeleton = connect_points_hmmf(endpoints, indices_fold,
            #                               likelihoods, neighbor_lists)
            if skeleton:

                #-------------------------------------------------------------
                # Smooth skeleton:
                #-------------------------------------------------------------
                if smooth_skeleton:
                    padded_skeleton = dilate(skeleton, nedges, neighbor_lists)
                    skeleton = connect_points_hmmf(endpoints, padded_skeleton,
                                                   likelihoods, neighbor_lists)

                fundi[skeleton] = fold_ID
                count += 1

    n_fundi = count
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.
          format(n_fundi, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, fundus points, and file name:
    #-------------------------------------------------------------------------
    fundi = fundi.tolist()

    if save_file:
        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(depth_file, fundi_file, fundi, 'fundi', folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundus_endpoints, fundi_file


# Example
if __name__ == "__main__" :

    # Extract fundus from one or more folds:
    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.plots import plot_vtk
    from mindboggle.features.fundi import extract_fundi
    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    single_fold = True
    if single_fold:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, True, True)
        fold_number = 11 #11
        folds[folds != fold_number] = -1
    else:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    #
    smooth_skeleton = False
    fundi, n_fundi, endpoints, fundi_file = extract_fundi(folds, depth_file,
        likelihoods_or_file, smooth_skeleton, save_file=True)
    #
    # View:
    plot_vtk(fundi_file)
