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
def extract_fundi(folds, sulci, curv_file, depth_file,
                  likelihoods, min_edges, min_distance,
                  erode_ratio, erode_min_size, smooth_skeleton, save_file):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a sulcus fold.

    Steps ::
        1. Find fundus endpoints (outer anchors) with find_outer_anchors().
        2. Include inner anchor points.
        3. Connect anchor points using connect_points_erosion().
        4. To do: Optionally smooth fundi.
        5. Segment fundi by sulcus definitions.

    Parameters
    ----------
    folds : list of integers
        fold number for each vertex
    curv_file :  string
        surface mesh file in VTK format with mean curvature values
    depth_file :  string
        surface mesh file in VTK format with rescaled depth values
    sulci : list of integers
        sulcus number for each vertex
    likelihoods : list of integers
        fundus likelihood value for each vertex
    min_edges : integer
        minimum number of edges between endpoints in find_outer_anchors()
    erosion_ratio : float
        fraction of indices to test for removal at each iteration
        in connect_points_erosion()
    smooth_skeleton : Boolean [Not yet implemented]
        smooth skeleton?
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    fundi : list of integers
        fundus numbers for all vertices (-1 for non-fundus vertices)
    n_fundi :  integer
        number of fundi
    fundi_file : string (if save_file)
        name of output VTK file with fundus numbers (-1 for non-fundus vertices)

    Examples
    --------
    >>> # Extract fundus from one or more folds:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.plots import plot_vtk
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(sulci_file, True, True)
    >>> curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'travel_depth_rescaled.vtk')
    >>> likelihoods_file = os.path.join(path, 'arno', 'shapes', 'likelihoods.vtk')
    >>> likelihoods, name = read_scalars(likelihoods_file, True, True)
    >>> single_fold = True
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> if single_fold:
    >>>     fold_number = 11 #11
    >>>     folds[folds != fold_number] = -1
    >>> min_edges = 10
    >>> min_distance = 10
    >>> erode_ratio = 0.10
    >>> erode_min_size = 10
    >>> smooth_skeleton = False
    >>> save_file = True
    >>> fundi, n_fundi, fundi_file = extract_fundi(folds, sulci, curv_file,
    >>>     depth_file, likelihoods, min_edges, min_distance,
    >>>     erode_ratio, erode_min_size, smooth_skeleton, save_file)
    >>> #
    >>> # View:
    >>> plot_vtk(fundi_file)

    """

    # Extract a skeleton to connect endpoints in a fold:
    import os
    import numpy as np
    from time import time

    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.compute import median_abs_dev
    from mindboggle.utils.paths import find_max_values
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.morph import dilate
    from mindboggle.utils.paths import find_outer_anchors, \
        connect_points_erosion, connect_points_hmmf

    # From connect_points_hmmf(): maximum neighborhood weight
    # (trust prior more for smoother fundi)
    wN_max = 2.0

    # Load values, threshold, and neighbors:
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
            # to serve as fundus endpoints :
            #-----------------------------------------------------------------
            outer_anchors, tracks = find_outer_anchors(indices_fold,
                neighbor_lists, values, depths, min_edges)

            #-----------------------------------------------------------------
            # Find inner anchor points:
            #-----------------------------------------------------------------
            inner_anchors = find_max_values(points, values, min_distance, thr)

            #-----------------------------------------------------------------
            # Connect endpoints to create skeleton:
            #-----------------------------------------------------------------
            B = -1 * np.ones(npoints)
            B[indices_fold] = 1
            skeleton = connect_points_erosion(B, neighbor_lists,
                outer_anchors, inner_anchors, values,
                erode_ratio, erode_min_size, save_steps=[], save_vtk='')
            if skeleton:

                #-------------------------------------------------------------
                # Smooth skeleton:
                #-------------------------------------------------------------
                if smooth_skeleton:

                    # Dilate the skeleton within the fold:
                    nedges = 2
                    print('    Dilate skeleton and intersect with fold...')
                    dilated = dilate(skeleton, nedges, neighbor_lists)
                    dilated = list(set(dilated).intersection(indices_fold))

                    # In case the dilation leads to topological changes,
                    # erode the fold again to the dilated skeleton (SLOW):
                    re_erode = False
                    if re_erode:
                        print('    Erode fold again to dilated skeleton...')
                        S = -1 * np.ones(npoints)
                        S[indices_fold] = 1
                        dilated = connect_points_erosion(S, neighbor_lists,
                            dilated, [], [], 1, 0, [], '')

                    # Set undilated likelihoods to -1 to preserve neighbors:
                    undilated = list(set(indices_fold).difference(dilated))
                    likelihoods_copy = likelihoods[:]
                    likelihoods_copy[undilated] = -1

                    # Smoothly re-skeletonize the dilated skeleton:
                    print('    Smoothly re-skeletonize dilated skeleton...')
                    skeleton = connect_points_hmmf(outer_anchors, dilated,
                        likelihoods_copy.tolist(), neighbor_lists, wN_max)

                    ## Plot overlap of dilated and pre-/post-smoothed skeleton:
                    #from mindboggle.utils.plots import plot_vtk
                    #D = -1*np.ones(npoints)
                    #D[dilated]=1; D[skeleton]=2; D[skeleton2]=3
                    #rewrite_scalars(curv_file, 'test.vtk', D, 'D', folds)
                    #plot_vtk('test.vtk')

                #-------------------------------------------------------------
                # Store skeleton:
                #-------------------------------------------------------------
                skeletons.extend(skeleton)

    #-------------------------------------------------------------------------
    # Create fundi by segmenting skeletons with overlapping sulcus labels:
    #-------------------------------------------------------------------------
    fundi = -1 * np.ones(npoints)
    indices = [x for x in skeletons if sulci[x] != -1]
    fundi[indices] = sulci[indices]

    n_fundi = len([x for x in np.unique(fundi) if x != -1])
    if n_fundi == 1:
        sdum = 'fundus'
    else:
        sdum = 'fundi'
    print('  ...Extracted {0} {1} ({2:.2f} seconds)'.
          format(n_fundi, sdum, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    fundi = fundi.tolist()

    if save_file:
        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(curv_file, fundi_file, fundi, 'fundi', folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundi_file


# Example
#if __name__ == "__main__" :

