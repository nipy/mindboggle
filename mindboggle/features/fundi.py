#!/usr/bin/env python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Arno Klein, 2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def extract_fundi(folds, curv_file, depth_file, min_separation=10,
                  erode_ratio=0.1, erode_min_size=1, save_file=False):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a fold.

    Steps ::
        1. Find fundus endpoints (outer anchors) with find_outer_anchors().
        2. Include inner anchor points.
        3. Connect anchor points using connect_points_erosion();
           inner anchors are removed if they result in endpoints.

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
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundus_per_fold : list of integers
        fundus numbers for all vertices, labeled by fold
        (-1 for non-fundus vertices)
    n_fundi_in_folds :  integer
        number of fundi
    fundus_per_fold_file : string (if save_file)
        output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> single_fold = True
    >>> import os
    >>> from mindboggle.io.vtks import read_scalars
    >>> from mindboggle.features.fundi import extract_fundi
    >>> from mindboggle.io.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
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
    >>> o1, o2, fundus_per_fold_file = extract_fundi(folds, curv_file,
    ...     depth_file, min_separation, erode_ratio, erode_min_size, save_file)
    >>> #
    >>> # View:
    >>> plot_surfaces(fundi_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np
    from time import time

    from mindboggle.io.vtks import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.guts.compute import median_abs_dev
    from mindboggle.guts.paths import find_max_values
    from mindboggle.guts.mesh import find_neighbors_from_file, find_complete_faces
    from mindboggle.guts.paths import find_outer_anchors, connect_points_erosion

    if isinstance(folds, list):
        folds = np.array(folds)

    # Load values, inner anchor threshold, and neighbors:
    faces, u1,u2, points, npoints, curvs, u3,u4 = read_vtk(curv_file, True,True)
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

            #-----------------------------------------------------------------
            # Remove fundus vertices if they complete triangle faces:
            #-----------------------------------------------------------------
            Iremove = find_complete_faces(skeletons, faces)
            if Iremove:
                skeletons = list(frozenset(skeletons).difference(Iremove))

    indices = [x for x in skeletons if folds[x] != -1]
    fundus_per_fold = -1 * np.ones(npoints)
    fundus_per_fold[indices] = folds[indices]
    n_fundi_in_folds = len([x for x in np.unique(fundus_per_fold)
                             if x != -1])
    if n_fundi_in_folds == 1:
        sdum = 'fold fundus'
    else:
        sdum = 'fold fundi'
    print('  ...Extracted {0} {1}; {2} total ({3:.2f} seconds)'.
          format(n_fundi_in_folds, sdum, n_fundi_in_folds, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    if n_fundi_in_folds > 0:
        fundus_per_fold = [int(x) for x in fundus_per_fold]
        if save_file:
            fundus_per_fold_file = os.path.join(os.getcwd(),
                                                'fundus_per_fold.vtk')
            rewrite_scalars(curv_file, fundus_per_fold_file, fundus_per_fold,
                            'fundi', folds)
            if not os.path.exists(fundus_per_fold_file):
                raise(IOError(fundus_per_fold_file + " not found"))
        else:
            fundus_per_fold_file = None

    return fundus_per_fold,  n_fundi_in_folds, fundus_per_fold_file


def segment_fundi(fundus_per_fold, sulci=[], vtk_file='', save_file=False):
    """
    Segment fundi by sulcus definitions.

    Parameters
    ----------
    fundus_per_fold : list of integers
        fundus numbers for all vertices, labeled by fold
        (-1 for non-fundus vertices)
    sulci : numpy array or list of integers
        sulcus number for each vertex, used to filter and label fundi
    vtk_file : string (if save_file)
        VTK file with sulcus number for each vertex
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundus_per_sulcus : list of integers
        fundus numbers for all vertices, labeled by sulcus
        (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundus_per_sulcus_file : string (if save_file)
        output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more sulci:
    >>> single_fold = True
    >>> import os
    >>> from mindboggle.io.vtks import read_scalars
    >>> from mindboggle.features.fundi import extract_fundi, segment_fundi
    >>> from mindboggle.io.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(vtk_file, True, True)
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
    >>> fundus_per_fold, o1, o2 = extract_fundi(folds, curv_file, depth_file, min_separation, erode_ratio, erode_min_size, save_file)
    >>> o1, o2, fundus_per_sulcus_file = segment_fundi(fundus_per_fold, sulci, vtk_file, save_file)
    >>> #
    >>> # View:
    >>> plot_surfaces(fundus_per_sulcus_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np

    from mindboggle.io.vtks import rewrite_scalars

    if isinstance(sulci, list):
        sulci = np.array(sulci)

    #-------------------------------------------------------------------------
    # Create fundi by segmenting fold fundi with overlapping sulcus labels:
    #-------------------------------------------------------------------------
    indices = [i for i,x in enumerate(fundus_per_fold) if x != -1]
    if indices and np.size(sulci):
        fundus_per_sulcus = -1 * np.ones(len(sulci))
        fundus_per_sulcus[indices] = sulci[indices]
        n_fundi = len([x for x in np.unique(fundus_per_sulcus) if x != -1])
    else:
        fundus_per_sulcus = []
        n_fundi = 0

    if n_fundi == 1:
        sdum = 'sulcus fundus'
    else:
        sdum = 'sulcus fundi'
    print('  Segmented {0} {1}'.format(n_fundi, sdum))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    fundus_per_sulcus_file = None
    if n_fundi > 0:
        fundus_per_sulcus = [int(x) for x in fundus_per_sulcus]
        if save_file and os.path.exists(vtk_file):
            fundus_per_sulcus_file = os.path.join(os.getcwd(),
                                                  'fundus_per_sulcus.vtk')
            rewrite_scalars(vtk_file, fundus_per_sulcus_file,
                            fundus_per_sulcus, 'fundus_per_sulcus',
                            fundus_per_sulcus)
            if not os.path.exists(fundus_per_sulcus_file):
                raise(IOError(fundus_per_sulcus_file + " not found"))

    return fundus_per_sulcus, n_fundi, fundus_per_sulcus_file
