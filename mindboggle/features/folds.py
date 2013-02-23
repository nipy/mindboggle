#!/usr/bin/env python
"""
Functions to extract folds.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#===============================================================================
# Extract folds
#===============================================================================
def extract_folds(depth_file, neighbor_lists=[], min_fold_size=1, extract_subfolds=True):
    """
    Use depth to extract folds from a triangular surface mesh.

    Steps ::
        1. Compute histogram of depth measures
        2. Find the deepest vertices
        3. Segment deep vertices as an initial set of folds
        4. Remove small folds
        5. Find and fill holes in the folds
        6. Optional: segment into subfolds
            a. Segment folds into "watershed basins"
            b. Shrink segments in folds with multiple segments
            c. Regrow shrunken segments
        7. Renumber segments

    Step 2::
        To extract an initial set of deep vertices from the surface mesh,
        we anticipate that there will be a rapidly decreasing distribution
        of low depth values (on the outer surface) with a long tail
        of higher depth values (in the folds), so we smooth the histogram's
        bin values (Gaussian), convolve to compute slopes,
        and find the depth value for the first bin with slope = 0.

    Step 5::
        The resulting separately numbered folds may have holes
        resulting from shallower areas within a fold, but calling fill_holes()
        at this stage can accidentally fill surfaces between folds,
        so we call fill_holes() with the argument exclude_values set to zero
        for zero-depth between folds.

    Step 6::
        The watershed() function has the same drawback as the segment() function
        -- the order of seed selection influences the result for multiple
        seeds within a connected set of vertices (region).
        To ameliorate this bias, we run the shrink_segments() function
        on the segments returned by watershed(), which shrinks segments in
        regions with multiple segments, and use these fractional segments
        as seeds for the propagate() function, which is slower and insensitive
        to depth, but is not biased by seed order.

    Parameters
    ----------
    depth_file : string
        surface mesh file in VTK format with faces and depth scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
        -- if empty list, construct from depth_file
    min_fold_size : int
        minimum fold size (number of vertices)
    extract_subfolds : Boolean
        segment folds into subfolds?

    Returns
    -------
    folds : array of integers
        fold numbers for all vertices (default -1 for non-fold vertices)
    n_folds :  int
        number of folds

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_faces_points, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.features.folds import extract_folds
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>>
    >>> folds, n_folds = extract_folds(depth_file, neighbor_lists, 50, True)
    >>>
    >>> # Write results to vtk file and view:
    >>> folds = folds.tolist()
    >>> rewrite_scalars(depth_file, 'test_extract_folds.vtk', folds, 'folds', folds)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_folds.vtk')

    """
    import numpy as np
    from time import time
    from scipy.ndimage.filters import gaussian_filter1d
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh import fill_holes
    from mindboggle.labels.segment import segment, watershed, shrink_segments
    if not len(neighbor_lists):
        from mindboggle.utils.mesh import find_neighbors

    do_fill_holes = True

    print("Extract folds in surface mesh")
    t0 = time()

    #---------------------------------------------------------------------------
    # Load depth values for all vertices
    #---------------------------------------------------------------------------
    faces, lines, indices, points, npoints, depths, \
        name = read_vtk(depth_file, return_first=True, return_array=True)

    #---------------------------------------------------------------------------
    # Find neighbors for each vertex
    #---------------------------------------------------------------------------
    if not neighbor_lists:
        neighbor_lists = find_neighbors(faces, npoints)

    #---------------------------------------------------------------------------
    # Compute histogram of depth measures
    #---------------------------------------------------------------------------
    min_vertices = 10000
    if npoints > min_vertices:
        nbins = np.round(npoints / 100.0)
    else:
        error("  Expecting at least {0} vertices to create depth histogram".
        format(min_vertices))
    bins, bin_edges = np.histogram(depths, bins=nbins)
    #>>> # Plot histogram:
    #>>> a,b,c = hist(depths, bins=nbins)

    #---------------------------------------------------------------------------
    # Anticipating that there will be a rapidly decreasing distribution
    # of low depth values (on the outer surface) with a long tail of higher
    # depth values (in the folds), smooth the bin values (Gaussian), convolve
    # to compute slopes, and find the depth for the first bin with slope = 0.
    #---------------------------------------------------------------------------
    bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    #>>> # Plot smoothed histogram:
    #>>> plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')
    window = [-1, 0, 1]
    bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
    ibins0 = np.where(bin_slopes == 0)[0]
    if ibins0.size:
        depth_threshold = bin_edges[ibins0[0]]
    else:
        depth_threshold = np.median(depths)

    #---------------------------------------------------------------------------
    # Find the deepest vertices
    #---------------------------------------------------------------------------
    indices_deep = [i for i,x in enumerate(depths) if x >= depth_threshold]
    if indices_deep:

        #-----------------------------------------------------------------------
        # Segment deep vertices as an initial set of folds
        #-----------------------------------------------------------------------
        print("  Segment vertices deeper than {0:.2f} as folds".format(depth_threshold))
        t1 = time()
        folds = segment(indices_deep, neighbor_lists)
        # Slightly slower alternative -- fill boundaries:
        #regions = -1 * np.ones(len(points))
        #regions[indices_deep] = 1
        #folds = segment_by_filling_boundaries(regions, neighbor_lists)
        print('    ...Segmented folds ({0:.2f} seconds)'.format(time() - t1))

        #-----------------------------------------------------------------------
        # Remove small folds
        #-----------------------------------------------------------------------
        if min_fold_size > 1:
            print('    Remove folds smaller than {0}'.format(min_fold_size))
            unique_folds = [x for x in np.unique(folds) if x > -1]
            for nfold in unique_folds:
                indices_fold = np.where(folds == nfold)[0]
                if len(indices_fold) < min_fold_size:
                    folds[indices_fold] = -1

        #-----------------------------------------------------------------------
        # Find and fill holes in the folds
        # Note: Surfaces surrounded by folds can be mistaken for holes,
        #       so exclude_values equals outer surface value of zero.
        #-----------------------------------------------------------------------
        if do_fill_holes:
            print("  Find and fill holes in the folds")
            folds = fill_holes(folds, neighbor_lists, exclude_values=[0],
                               values=depths)

        #-----------------------------------------------------------------------
        # Extract subfolds
        #-----------------------------------------------------------------------
        if extract_subfolds:

            #-------------------------------------------------------------------
            # Segment folds into "watershed basins"
            #-------------------------------------------------------------------
            indices_folds = [i for i,x in enumerate(folds) if x > -1]
            #indices_folds = np.where(folds > -1)[0]
            segments = watershed(depths, indices_folds, neighbor_lists,
                                 depth_ratio=0.1, tolerance=0.01)

            #-------------------------------------------------------------------
            # Shrink segments in folds with multiple segments
            #-------------------------------------------------------------------
            shrunken_segments = shrink_segments(folds, segments, depths,
                                                remove_fraction=0.75,
                                                only_multiple_segments=True)
            print('  ...Segmented, removed small folds, filled holes, shrank segments'
                  ' ({0:.2f} seconds)'.format(time() - t0))

            #-------------------------------------------------------------------
            # Regrow shrunken segments
            #-------------------------------------------------------------------
            print("  Regrow shrunken segments")
            t2 = time()
            seed_lists = []
            unique_shrunken = [x for x in np.unique(shrunken_segments) if x > -1]
            for n_shrunken in unique_shrunken:
                seed_lists.append([i for i,x in enumerate(shrunken_segments)
                                   if x == n_shrunken])
            regrown_segments = segment(indices_folds, neighbor_lists,
                                       min_fold_size, seed_lists)
            #regrown_segments = propagate(points, faces, folds, shrunken_segments,
            #                             folds, max_iters=1000, tol=0.001, sigma=5)
            folds[regrown_segments > -1] = regrown_segments[regrown_segments > -1] + \
                                           np.max(folds) + 1
            print('    ...Segmented individual folds ({0:.2f} seconds)'.
                  format(time() - t2))

        #-----------------------------------------------------------------------
        # Renumber folds so they are sequential
        #-----------------------------------------------------------------------
        renumber_folds = -1 * np.ones(len(folds))
        fold_numbers = [int(x) for x in np.unique(folds) if x > -1]
        for i_fold, n_fold in enumerate(fold_numbers):
            fold = [i for i,x in enumerate(folds) if x == n_fold]
            renumber_folds[fold] = i_fold
        folds = renumber_folds
        n_folds = i_fold + 1

        # Print statement
        print('  ...Extracted {0} folds ({1:.2f} seconds)'.
              format(n_folds, time() - t0))
    else:
        print('  No deep vertices')

    #---------------------------------------------------------------------------
    # Return folds, number of folds
    #---------------------------------------------------------------------------
    return folds, n_folds


#===============================================================================
# Normalize depths in each fold
#===============================================================================
def normalize_fold_depths(depth_file, folds, save_file=False):

    """
    Normalize depths in each fold.

    Parameters
    ----------
    depth_file : string
        name of VTK file with a depth value for each vertex
    folds : list
        fold ID for each vertex

    Returns
    -------
    norm_depth_folds : list of floats
        depth values for fold IDs >-1, normalized for each fold
    save_file : Boolean
        save output VTK file?
    norm_depth_file : string (if save_file)

        name of output VTK file with normalized depth values in each fold

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.features.folds import normalize_fold_depths
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> #
    >>> normalize_fold_depths(depth_file, folds, save_file=True)
    >>> # View:
    >>> norm_depth_file = os.path.join(os.getcwd(),
    >>>     os.path.basename(depth_file).strip('vtk') + 'norm.vtk')
    >>> os.system('mayavi2 -m Surface -d ' + norm_depth_file + '&')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    depths, name = read_scalars(depth_file, True, True)


    norm_depth_folds = -1 * np.ones(len(folds))

    unique_fold_IDs = np.unique(folds)
    unique_fold_IDs = [x for x in unique_fold_IDs if x >= 0]

    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:

            # Normalize fold depth values
            norm_depth_folds[indices_fold] = depths[indices_fold] / np.max(depths[indices_fold])

    if save_file:

        norm_depth_file = os.path.join(os.getcwd(),
                                       os.path.basename(depth_file).strip('vtk') + 'norm.vtk')
        rewrite_scalars(depth_file, norm_depth_file,
                        norm_depth_folds, 'norm_depths', norm_depth_folds)

        return norm_depth_folds, norm_depth_file
    else:
        return norm_depth_folds
