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
def extract_fundi(folds, sulci, likelihoods, rescaled_depth_file,
                  depth_file, min_edges=10, erosion_ratio=0.25,
                  smooth_skeleton=False, filter=False, filter_file='',
                  save_file=False):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a sulcus fold.

    Steps ::

        1. Find fundus endpoints using find_outer_anchors().
        2. Filter out large-value vertices by reducing their likelihoods.
        3. Connect fundus endpoints using connect_points_erosion().
        4. To do: Optionally smooth fundi.
        5. Segment fundi by sulcus definitions.

    Parameters
    ----------
    folds : list of integers
        fold number for each vertex
    sulci : list of integers
        sulcus number for each vertex
    likelihoods : list of integers
        fundus likelihood value for each vertex
    rescaled_depth_file :  string
        surface mesh file in VTK format with scalar rescaled depth values
    depth_file :  string
        surface mesh file in VTK format with (complete) scalar depth values
    min_edges : integer
        minimum number of edges between endpoints in find_outer_anchors()
    erosion_ratio : float
        fraction of indices to test for removal at each iteration
        in connect_points_erosion()
    smooth_skeleton : Boolean [Not yet implemented]
        smooth skeleton?
    filter : Boolean
        filter vertices to connect based on values in a file?
    filter_file : list of integers
        values used to filter vertices by reducing their likelihoods
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
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.plots import plot_vtk
    >>> from mindboggle.features.fundi import extract_fundi
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> sulci, name = read_scalars(sulci_file, True, True)
    >>> likelihoods_file = os.path.join(path, 'arno', 'shapes', 'likelihoods.vtk')
    >>> likelihoods, name = read_scalars(likelihoods_file, True, True)
    >>> rescaled_depth_file = os.path.join(path, 'arno', 'shapes', 'travel_depth_rescaled.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> single_fold = False
    >>> if single_fold:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, True, True)
    >>>     fold_number = 11 #11
    >>>     folds[folds != fold_number] = -1
    >>> else:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, True, True)
    >>> #
    >>> min_edges = 10
    >>> erosion_ratio = 0.25
    >>> smooth_skeleton = True
    >>> filter = False
    >>> filter_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> save_file = True
    >>> fundi, n_fundi, fundi_file = extract_fundi(folds, sulci, likelihoods,
    >>>     rescaled_depth_file, depth_file, min_edges, erosion_ratio,
    >>>     smooth_skeleton, filter, filter_file, save_file)
    >>> #
    >>> # View:
    >>> plot_vtk(fundi_file)

    """
    import os
    import numpy as np
    from time import time

    from mindboggle.utils.paths import find_outer_anchors, \
        connect_points_erosion, connect_points_hmmf
    from mindboggle.utils.mesh import find_neighbors_from_file, find_neighbors
    from mindboggle.utils.morph import dilate
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Run connect_points_erosion() or connect_points_hmmf():
    run_erosion = True

    # From connect_points_hmmf():
    # maximum neighborhood weight (trust prior more for smoother fundi)
    wN_max = 2.0

    # Load depths and neighbors:
    neighbor_lists = find_neighbors_from_file(depth_file)
    depths, name = read_scalars(rescaled_depth_file)
    npoints = len(depths)

    #-------------------------------------------------------------------------
    # Loop through folds:
    #-------------------------------------------------------------------------
    t1 = time()
    skeletons = []
    unique_fold_IDs = [x for x in np.unique(folds) if x != -1]

    print("Extract a fundus from each of {0} regions...".
          format(len(unique_fold_IDs)))

    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:
            print('  Fold {0}:'.format(int(fold_ID)))

            #-----------------------------------------------------------------
            # Find fundus endpoints on the boundary of the surface region:
            #-----------------------------------------------------------------
            endpoints, endtracks = find_outer_anchors(indices_fold,
                neighbor_lists, likelihoods, depths, min_edges)

            #-----------------------------------------------------------------
            # Filter out large-value vertices by reducing their likelihoods:
            #-----------------------------------------------------------------
            if filter:
                values, name = read_scalars(filter_file, True, True)
                V = [x for x in values if x != -1]
                threshold = np.median(V) + np.std(V)
                likelihoods[np.where(np.abs(values) > threshold)[0]] /= 10

            #-----------------------------------------------------------------
            # Connect endpoints to create skeleton:
            #-----------------------------------------------------------------
            if run_erosion:
                B = -1 * np.ones(npoints)
                B[indices_fold] = 1
                skeleton = connect_points_erosion(B, endpoints, neighbor_lists,
                                                  likelihoods, erosion_ratio)
            else:
                skeleton = connect_points_hmmf(endpoints, indices_fold,
                                               likelihoods, neighbor_lists,
                                               wN_max=2.0)
            if skeleton:

                #-------------------------------------------------------------
                # Smooth skeleton:
                #-------------------------------------------------------------
                if run_erosion and smooth_skeleton:

                    # Dilate the skeleton within the fold:
                    nedges = 2
                    print('    Dilate skeleton and intersect with fold...')
                    dilated = dilate(skeleton, nedges, neighbor_lists)
                    dilated = list(frozenset(dilated).intersection(indices_fold))

                    # In case the dilation leads to topological changes,
                    # erode the fold again to the dilated skeleton (SLOW):
                    re_erode = False
                    if re_erode:
                        print('    Erode fold again to dilated skeleton...')
                        S = -1 * np.ones(npoints)
                        S[indices_fold] = 1
                        dilated = connect_points_erosion(S, dilated,
                            neighbor_lists, values=[], erosion_ratio=1)

                    # Set undilated likelihoods to -1 to preserve neighbors:
                    undilated = list(frozenset(indices_fold).difference(dilated))
                    likelihoods_copy = likelihoods[:]
                    likelihoods_copy[undilated] = -1

                    # Smoothly re-skeletonize the dilated skeleton:
                    print('    Smoothly re-skeletonize dilated skeleton...')
                    skeleton = connect_points_hmmf(endpoints, dilated,
                        likelihoods_copy.tolist(), neighbor_lists, wN_max)

                    # Plot overlap of dilated and pre-/post-smoothed skeleton:
                    #skeleton2 = connect_points_hmmf(endpoints, dilated,
                    #    likelihoods_copy, neighbor_lists)
                    #D = -1*np.ones(npoints)
                    #D[dilated]=1; D[skeleton]=2; D[skeleton2]=3
                    #rewrite_scalars(depth_file, 'test.vtk', D, 'D', folds)
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
    print('  ...Extracted {0} fundi ({1:.2f} seconds)'.
          format(n_fundi, time() - t1))

    #-------------------------------------------------------------------------
    # Return fundi, number of fundi, and file name:
    #-------------------------------------------------------------------------
    fundi = fundi.tolist()

    if save_file:
        fundi_file = os.path.join(os.getcwd(), 'fundi.vtk')
        rewrite_scalars(depth_file, fundi_file, fundi, 'fundi', folds)
    else:
        fundi_file = None

    return fundi, n_fundi, fundi_file


# Example
if __name__ == "__main__" :

    # Extract fundus from one or more folds:
    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.plots import plot_vtk
    from mindboggle.features.fundi import extract_fundi
    path = os.environ['MINDBOGGLE_DATA']
    sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    sulci, name = read_scalars(sulci_file, True, True)
    likelihoods_file = os.path.join(path, 'arno', 'shapes', 'likelihoods.vtk')
    likelihoods, name = read_scalars(likelihoods_file, True, True)
    rescaled_depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    single_fold = True
    if single_fold:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, True, True)
        fold_number = 11 #11
        folds[folds != fold_number] = -1
    else:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, True, True)
    #
    min_edges = 10
    erosion_ratio = 0.25
    smooth_skeleton = True
    filter = False
    filter_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    save_file = True
    fundi, n_fundi, fundi_file = extract_fundi(folds, sulci, likelihoods,
        rescaled_depth_file, depth_file, min_edges, erosion_ratio,
        smooth_skeleton, filter, filter_file, save_file)
    #
    # View:
    plot_vtk(fundi_file)
