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
def extract_fundi(folds, curv_file, depth_file, min_separation=10,
                  erode_ratio=0.1, erode_min_size=1, sulci=[],
                  save_file=False):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a sulcus fold.

    Steps ::
        1. Find fundus endpoints (outer anchors) with find_outer_anchors().
        2. Include inner anchor points.
        3. Connect anchor points using connect_points_erosion();
           inner anchors are removed if they result in endpoints.
        4. Optionally segment fundi by sulcus definitions.
           Return both fundi in all folds as well as those in sulcus folds.

        Optional postprocessing step: smooth with smooth_skeleton().

    Parameters
    ----------
    folds : numpy array or list of integers
        fold number for each vertex
    curv_file :  string
        surface mesh file in VTK format with mean curvature values
    depth_file :  string
        surface mesh file in VTK format with rescaled depth values
    likelihoods : list of integers
        fundus likelihood value for each vertex
    min_separation : integer
        minimum number of edges between inner/outer anchor points
    erode_ratio : float
        fraction of indices to test for removal at each iteration
        in connect_points_erosion()
    sulci : numpy array or list of integers
        sulcus number for each vertex, used to filter and label fundi
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundi : list of integers
        fundus numbers for all vertices, labeled by sulcus
        (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundi_file : string (if save_file)
        output VTK file with fundus numbers (-1 for non-fundus vertices)
    fundi_all_folds : list of integers
        same as fundi, but labeled by all folds rather than just by sulci
    n_fundi_all_folds :  integer
        same as n_fundi, but for all folds
    fundi_all_folds_file : string (if save_file)
        same as fundi_file, but for all folds

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> single_fold = True
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.features.fundi import extract_fundi
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(sulci_file, True, True)
    >>> curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'travel_depth_rescaled.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> if single_fold:
    >>>     fold_number = 2 #11
    >>>     folds[folds != fold_number] = -1
    >>> min_separation = 10
    >>> erode_ratio = 0.10
    >>> erode_min_size = 10
    >>> save_file = True
    >>> o1, o2, fundi_file, o3, o4, o5 = extract_fundi(folds, curv_file,
    ...     depth_file, min_separation, erode_ratio, erode_min_size, sulci,
    ...     save_file)
    >>> #
    >>> # View:
    >>> plot_surfaces(fundi_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np
    from time import time

    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.compute import median_abs_dev
    from mindboggle.utils.paths import find_max_values
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.paths import find_outer_anchors, connect_points_erosion

    if isinstance(folds, list):
        folds = np.array(folds)
    if isinstance(sulci, list):
        sulci = np.array(sulci)

    # Load values, inner anchor threshold, and neighbors:
    u1,u2,u3, points, npoints, curvs, u4,u5 = read_vtk(curv_file, True,True)
    depths, name = read_scalars(depth_file, True, True)
    values = curvs * depths
    values0 = [x for x in values if x > 0]
    thr = np.median(values0) + 2 * median_abs_dev(values0)
    neighbor_lists = find_neighbors_from_file(curv_file)

    #-------------------------------------------------------------------------
    # Loop through folds:
    #-------------------------------------------------------------------------
    t1 = time()
    skeletons = []
    unique_fold_IDs = [x for x in np.unique(folds) if x != -1]

    if len(unique_fold_IDs) == 1:
        print("Extract a fundus from 1 fold...")
    else:
        print("Extract a fundus from each of {0} folds...".
              format(len(unique_fold_IDs)))

    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Fold {0}:'.format(int(fold_ID)))

            #-----------------------------------------------------------------
            # Find outer anchor points on the boundary of the surface region,
            # to serve as fundus endpoints:
            #-----------------------------------------------------------------
            outer_anchors, tracks = find_outer_anchors(indices_fold,
                neighbor_lists, values, depths, min_separation)

            #-----------------------------------------------------------------
            # Find inner anchor points:
            #-----------------------------------------------------------------
            inner_anchors = find_max_values(points, values, min_separation, thr)

            #-----------------------------------------------------------------
            # Connect anchor points to create skeleton:
            #-----------------------------------------------------------------
            B = -1 * np.ones(npoints)
            B[indices_fold] = 1
            skeleton = connect_points_erosion(B, neighbor_lists,
                outer_anchors, inner_anchors, values,
                erode_ratio, erode_min_size, save_steps=[], save_vtk='')
            if skeleton:
                skeletons.extend(skeleton)

    #-------------------------------------------------------------------------
    # Create fundi by segmenting skeletons with overlapping sulcus labels:
    #-------------------------------------------------------------------------
    indices = [x for x in skeletons if folds[x] != -1]
    fundi_all_folds = -1 * np.ones(npoints)
    fundi_all_folds[indices] = folds[indices]
    n_fundi_all_folds = len([x for x in np.unique(fundi_all_folds)
                             if x != -1])
    if np.size(sulci):
        indices = [x for x in skeletons if sulci[x] != -1]
        fundi = -1 * np.ones(npoints)
        fundi[indices] = sulci[indices]
        n_fundi = len([x for x in np.unique(fundi) if x != -1])
    else:
        fundi = []
        fundi_file = None
        n_fundi = 0

    if n_fundi == 1:
        sdum = 'fundus'
    else:
        sdum = 'fundi'
    print('  ...Extracted {0} {1}; {2} total ({3:.2f} seconds)'.
          format(n_fundi, sdum, n_fundi_all_folds, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    if n_fundi > 0:
        fundi = fundi.tolist()
        if save_file:
            fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
            rewrite_scalars(curv_file, fundi_file, fundi, 'fundi', folds)
            if not os.path.exists(fundi_file):
                raise(IOError(fundi_file + " not found"))
        else:
            fundi_file = None

    if n_fundi_all_folds > 0:
        fundi_all_folds = fundi_all_folds.tolist()
        if save_file:
            fundi_all_folds_file = os.path.join(os.getcwd(),
                                                'fundi_all_folds.vtk')
            rewrite_scalars(curv_file, fundi_all_folds_file, fundi_all_folds,
                            'fundi', folds)
            if not os.path.exists(fundi_all_folds_file):
                raise(IOError(fundi_all_folds_file + " not found"))
        else:
            fundi_all_folds_file = None

    return fundi, n_fundi, fundi_file, \
           fundi_all_folds,  n_fundi_all_folds, fundi_all_folds_file

