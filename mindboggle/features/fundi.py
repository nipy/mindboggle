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
                  depth_file, smooth_skeleton=False, save_file=False):
    """
    Extract fundi from folds.

    A fundus is a branching curve that runs along the deepest and most
    highly curved portions of a sulcus fold.

    Steps ::

        1. Find fundus endpoints using find_outer_anchors().
        2. Connect fundus endpoints using connect_points_erosion().
        3. FIX: Optionally smooth fundi using connect_points_hmmf().
        4. Segment fundi by sulcus definitions.

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
    >>> likelihoods_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> likelihoods, name = read_scalars(likelihoods_file, True, True)
    >>> rescaled_depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> single_fold = True
    >>> if single_fold:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, True, True)
    >>>     fold_number = 11 #11
    >>>     folds[folds != fold_number] = -1
    >>> else:
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     folds, name = read_scalars(folds_file, True, True)
    >>> #
    >>> smooth_skeleton = False
    >>> save_file = True
    >>> fundi, n_fundi, fundi_file = extract_fundi(folds, sulci, likelihoods,
    >>>     rescaled_depth_file, depth_file, smooth_skeleton, save_file)
    >>> #
    >>> # View:
    >>> plot_vtk(fundi_file)

    """
    import os
    import numpy as np
    from time import time

    from mindboggle.utils.paths import find_outer_anchors, \
        connect_points_erosion, connect_points_hmmf
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.utils.morph import dilate
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Load depths and neighbors:
    neighbor_lists = find_neighbors_from_file(depth_file)
    depths, name = read_scalars(rescaled_depth_file)

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
                neighbor_lists, likelihoods, depths, min_edges=5)

            #-----------------------------------------------------------------
            # Connect endpoints to create skeleton:
            #-----------------------------------------------------------------
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
                    nedges = 2
                    padded_skeleton = dilate(skeleton, nedges, neighbor_lists)
                    skeleton = connect_points_hmmf(endpoints, padded_skeleton,
                                                   likelihoods, neighbor_lists)

                skeletons.extend(skeleton)

    #-------------------------------------------------------------------------
    # Create fundi by segmenting skeletons with overlapping sulcus labels:
    #-------------------------------------------------------------------------
    fundi = -1 * np.ones(len(folds))
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
    likelihoods_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    likelihoods, name = read_scalars(likelihoods_file, True, True)
    rescaled_depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    single_fold = True
    if single_fold:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, True, True)
        fold_number = 1 #11
        folds[folds != fold_number] = -1
    else:
        folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
        folds, name = read_scalars(folds_file, True, True)
    #
    smooth_skeleton = False
    save_file = True
    fundi, n_fundi, fundi_file = extract_fundi(folds, sulci, likelihoods,
        rescaled_depth_file, depth_file, smooth_skeleton, save_file)
    #
    # View:
    plot_vtk(fundi_file)
