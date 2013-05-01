#!/usr/bin/env python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame, 2012-2013  .  yrjo.hame@gmail.com
Arno Klein, 2012-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Extract all fundi
#=============================================================================
def extract_fundi(folds_or_file, depth_file, likelihoods_or_file,
                  smooth_skeleton=False, save_file=False):
    """
    Extract fundi from folds.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Steps ::

        1. Find fundus endpoints.
        2. Connect fundus endpoints and extract fundi.

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
    from mindboggle.utils.mesh import find_neighbors_from_file, dilate
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

    # Extract fundus from a single fold or multiple folds:

    """
    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.plots import plot_vtk
    from mindboggle.features.fundi import extract_fundi

    path = os.environ['MINDBOGGLE_DATA']
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    likelihoods_or_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')

    single_fold = False
    if single_fold:
        folds_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    else:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    folds, name = read_scalars(folds_file, return_first=True, return_array=True)

    fundi, n_fundi, endpoints, fundi_file = extract_fundi(folds, depth_file,
        likelihoods_or_file, save_file=True)

    # View:
    plot_vtk(fundi_file)
    """

    # Setup:
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.plots import plot_vtk
    from mindboggle.utils.paths import track_endpoints

    single_fold = False

    path = os.environ['MINDBOGGLE_DATA']
    values_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    values_seeding_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    values, name = read_scalars(values_file, True, True)
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    values_seeding, name = read_scalars(values_seeding_file, True, True)
    neighbor_lists = find_neighbors_from_file(depth_file)
    min_size = 50
    min_edges = 5
    use_threshold = True

    #-------------------------------------------------------------------------
    # Extract tracks and endpoints on a single fold:
    #-------------------------------------------------------------------------
    if single_fold:
        fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
        fold, name = read_scalars(fold_file)
        indices = [i for i,x in enumerate(fold) if x != -1]

        tracks, indices_endpoints = track_endpoints(indices, neighbor_lists, \
            values, values_seeding, min_edges, use_threshold, backtrack=True)

        # View results atop values:
        indices_tracks = [x for lst in tracks for x in lst]
        values[indices_tracks] = max(values) + 0.5
        values[indices_endpoints] = max(values) + 0.1
        rewrite_scalars(depth_file, 'track_endpoints.vtk', \
                        values, 'endpoints_on_values_in_fold', fold)
        plot_vtk('track_endpoints.vtk')
    #-------------------------------------------------------------------------
    # Extract tracks and endpoints on every fold in a hemisphere:
    #-------------------------------------------------------------------------
    else:
        plot_each_fold = False
        folds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
        folds, name = read_scalars(folds_file)
        fold_numbers = [x for x in np.unique(folds) if x != -1]
        nfolds = len(fold_numbers)
        all_endpoints = []
        all_tracks = []
        for ifold, fold_number in enumerate(fold_numbers):
            print('Fold {0} ({1} of {2})'.format(int(fold_number), ifold+1, nfolds))
            indices = [i for i,x in enumerate(folds) if x == fold_number]
            if len(indices) > min_size:
                tracks, indices_endpoints = track_endpoints(indices,
                    neighbor_lists, values, values_seeding, min_edges,
                    use_threshold, backtrack=True)
                if tracks:
                    indices_tracks = [x for lst in tracks for x in lst]
                    all_endpoints.extend(indices_endpoints)
                    all_tracks.extend(indices_tracks)
                    # Plot each fold:
                    if plot_each_fold:
                        fold = -1 * np.ones(len(values))
                        fold[indices] = 1
                        values[indices_tracks] = max(values) + 0.5
                        values[indices_endpoints] = max(values) + 0.1
                        rewrite_scalars(depth_file, 'track_endpoints.vtk',
                                values, 'endpoints_on_values_in_fold', fold)
                        plot_vtk('track_endpoints.vtk')
        T = -1 * np.ones(len(values))
        T[all_tracks] = 1
        T[all_endpoints] = 2

        # Write results to VTK file and view:
        rewrite_scalars(folds_file, 'track_endpoints.vtk',
                        T, 'tracks_endpoints_on_folds', folds)
        plot_vtk('track_endpoints.vtk')
