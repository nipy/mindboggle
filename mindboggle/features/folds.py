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
def extract_folds(depth_file, min_fold_size=50, tiny_depth=0.001, save_file=False):
    """
    Use depth to extract folds from a triangular surface mesh.

    Steps ::
        1. Compute histogram of depth measures.
        2. Define a depth threshold and find the deepest vertices.
        3. Segment deep vertices as an initial set of folds.
        4. Remove small folds.
        5. Find and fill holes in the folds.
        6. Renumber folds.

    Step 2 ::
        To extract an initial set of deep vertices from the surface mesh,
        we anticipate that there will be a rapidly decreasing distribution
        of low depth values (on the outer surface) with a long tail
        of higher depth values (in the folds), so we smooth the histogram's
        bin values, convolve to compute slopes, and find the depth value
        for the first bin with slope = 0. This is our threshold.

    Step 5 ::
        The folds could have holes in areas shallower than the depth threshold.
        Calling fill_holes() could accidentally include very shallow areas
        (in an annulus-shaped fold, for example), so we call fill_holes() with
        the argument exclude_range set close to zero to retain these areas.

    Parameters
    ----------
    depth_file : string
        surface mesh file in VTK format with faces and depth scalar values
    min_fold_size : int
        minimum fold size (number of vertices)
    tiny_depth : float
        largest non-zero depth value that will stop a hole from being filled
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    folds : list of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    n_folds :  int
        number of folds
    depth_threshold :  float
        threshold defining the minimum depth for vertices to be in a fold
    bins :  list of integers
        histogram bins: each is the number of vertices within a range of depth values
    bin_edges :  list of floats
        histogram bin edge values defining the bin ranges of depth values
    folds_file : string (if save_file)
        name of output VTK file with fold IDs (-1 for non-fold vertices)

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> import pylab
    >>> from scipy.ndimage.filters import gaussian_filter1d
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file, plot_vtk
    >>> from mindboggle.features.folds import extract_folds
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> min_fold_size = 50
    >>> tiny_depth = 0.001
    >>> save_file = True
    >>> folds, n_folds, thr, bins, bin_edges, folds_file = extract_folds(depth_file,
    >>>     min_fold_size, tiny_depth, save_file)
    >>> # View folds:
    >>> plot_vtk('folds.vtk')
    >>> # Plot histogram and depth threshold:
    >>> depths, name = read_scalars(depth_file)
    >>> nbins = np.round(len(depths) / 100.0)
    >>> a,b,c = pylab.hist(depths, bins=nbins)
    >>> pylab.plot(thr*np.ones((100,1)),
    >>>            np.linspace(0, max(bins), 100), 'r.')
    >>> pylab.show()
    >>> # Plot smoothed histogram:
    >>> bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    >>> pylab.plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')
    >>> pylab.show()

    """
    import os
    import sys
    import numpy as np
    from time import time
    from scipy.ndimage.filters import gaussian_filter1d
    from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk
    from mindboggle.utils.mesh import find_neighbors, fill_holes
    from mindboggle.labels.segment import segment

    do_fill_holes = True

    print("Extract folds in surface mesh")
    t0 = time()

    #---------------------------------------------------------------------------
    # Load depth values for all vertices
    #---------------------------------------------------------------------------
    faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
        return_first=True, return_array=True)

    #---------------------------------------------------------------------------
    # Find neighbors for each vertex
    #---------------------------------------------------------------------------
    neighbor_lists = find_neighbors(faces, npoints)

    #---------------------------------------------------------------------------
    # Compute histogram of depth measures
    #---------------------------------------------------------------------------
    min_vertices = 10000
    if npoints > min_vertices:
        nbins = np.round(npoints / 100.0)
    else:
        sys.err("  Expecting at least {0} vertices to create depth histogram".
            format(min_vertices))
    bins, bin_edges = np.histogram(depths, bins=nbins)

    #---------------------------------------------------------------------------
    # Anticipating that there will be a rapidly decreasing distribution
    # of low depth values (on the outer surface) with a long tail of higher
    # depth values (in the folds), smooth the bin values (Gaussian), convolve
    # to compute slopes, and find the depth for the first bin with slope = 0.
    #---------------------------------------------------------------------------
    bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    window = [-1, 0, 1]
    bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
    ibins0 = np.where(bin_slopes == 0)[0]
    if ibins0.shape:
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
        print('  ...Segmented folds ({0:.2f} seconds)'.format(time() - t1))

        #-----------------------------------------------------------------------
        # Remove small folds
        #-----------------------------------------------------------------------
        if min_fold_size > 1:
            print('  Remove folds smaller than {0}'.format(min_fold_size))
            unique_folds = [x for x in np.unique(folds) if x > -1]
            for nfold in unique_folds:
                indices_fold = np.where(folds == nfold)[0]
                if len(indices_fold) < min_fold_size:
                    folds[indices_fold] = -1

        #-----------------------------------------------------------------------
        # Find and fill holes in the folds
        # Note: Surfaces surrounded by folds can be mistaken for holes,
        #       so exclude_range includes outer surface values close to zero.
        #-----------------------------------------------------------------------
        if do_fill_holes:
            print("  Find and fill holes in the folds")
            folds = fill_holes(folds, neighbor_lists, values=depths,
                               exclude_range=[0, tiny_depth])

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
    # Return folds, number of folds, file name
    #---------------------------------------------------------------------------
    if save_file:

        folds_file = os.path.join(os.getcwd(), 'folds.vtk')
        rewrite_scalars(depth_file, folds_file, folds, 'folds', folds)

    else:
        folds_file = None

    return folds.tolist(), n_folds, depth_threshold, bins, bin_edges, folds_file


#===============================================================================
# Extract subfolds
#===============================================================================
def extract_subfolds(depth_file, folds, min_subfold_size=50, depth_factor=0.25,
                     depth_ratio=0.1, tolerance=0.01, shrink_factor=0.5,
                     save_file=False):
    """
    Use depth to segment folds into subfolds in a triangular surface mesh.

    Steps ::
        1. Segment folds into "watershed basins".
        2. Shrink segments by either:
            a. removing shallow vertices (in folds with multiple segments), or
            b. using the watershed basin seeds
        3. Regrow shrunken segments as subfolds.
        4. Renumber subfolds.

    Step 1 ::
        The watershed() function, when used alone, has the same drawback as
        segment() -- the order of seed selection influences the result for
        multiple seeds within a connected set of vertices (region).
        To ameliorate this bias, we run the shrink_segments() function
        on the segments returned by watershed(), which shrinks segments in
        regions with multiple segments, and use these fractional segments
        as seeds for the segment() function, or the propagate() function,
        which is slower and insensitive to depth, but is not biased by seed order.

    Parameters
    ----------
    depth_file : string
        surface mesh file in VTK format with faces and depth scalar values
    folds : list of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    min_subfold_size : int
        minimum fold size (number of vertices)
    depth_factor : float
        watershed() depth_factor:
        factor to determine whether to merge two neighboring watershed catchment
        basins -- they are merged if the Euclidean distance between their basin
        seeds is less than this fraction of the maximum Euclidean distance
        between points having minimum and maximum depths
    depth_ratio : float
        watershed() depth_ratio:
        the minimum fraction of depth for a neighboring shallower
        watershed catchment basin (otherwise merged with the deeper basin)
    tolerance : float
        watershed() tolerance:
        tolerance for detecting differences in depth between vertices
    shrink_factor : float
        shrink each fold to this fraction of its maximum depth
        for folds with multiple subfolds, to regrow the subfolds
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    subfolds : list of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    n_subfolds :  int
        number of subfolds
    subfolds_file : string (if save_file)
        name of output VTK file with fold IDs (-1 for non-fold vertices)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.folds import extract_subfolds
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> min_subfold_size = 50
    >>> depth_factor = 0.5
    >>> depth_ratio = 0.1
    >>> tolerance = 0.01
    >>> shrink_factor = 0.5
    >>> #
    >>> subfolds, n_subfolds, subfolds_file = extract_subfolds(depth_file,
    >>>     folds, min_subfold_size, depth_factor, depth_ratio, tolerance,
    >>>     shrink_factor, True)
    >>> #
    >>> # View:
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('subfolds.vtk')

    """
    import os
    import numpy as np
    from time import time
    from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.labels.segment import segment, propagate, watershed, \
                                          shrink_segments

    # Shrink folds or use seeds:
    regrow = True
    regrow_shrunken = True

    # Use propagate() to regrow subfolds:
    regrow_by_propagating = 0#True

    print("Segment folds into subfolds")
    t0 = time()

    #---------------------------------------------------------------------------
    # Load depth values for all vertices
    #---------------------------------------------------------------------------
    faces, lines, indices, points, npoints, depths, \
        name = read_vtk(depth_file, return_first=True, return_array=True)

    #---------------------------------------------------------------------------
    # Find neighbors for each vertex
    #---------------------------------------------------------------------------
    neighbor_lists = find_neighbors(faces, npoints)

    #---------------------------------------------------------------------------
    # Segment folds into "watershed basins"
    #---------------------------------------------------------------------------
    indices_folds = [i for i,x in enumerate(folds) if x > -1]
    subfolds, seed_indices = watershed(depths, points, indices_folds,
                                       neighbor_lists, depth_factor,
                                       depth_ratio, tolerance)

    #---------------------------------------------------------------------------
    # Regrow segments
    #---------------------------------------------------------------------------
    if regrow:

        #-----------------------------------------------------------------------
        # Regrow segments after shrinking them in depth
        #-----------------------------------------------------------------------
        if regrow_shrunken:
            print("  Regrow segments after shrinking in depth")
            shrunken = shrink_segments(folds, subfolds, depths,
                                       shrink_factor=shrink_factor,
                                       only_multiple_segments=True)
            if not regrow_by_propagating:
                seed_lists = []
                unique_shrunken = [x for x in np.unique(shrunken)
                                   if x > -1]
                for n_shrunken in unique_shrunken:
                    seed_lists.append([i for i,x in enumerate(shrunken)
                                       if x == n_shrunken])
        #-----------------------------------------------------------------------
        # Regrow segments from watershed basin seeds
        #-----------------------------------------------------------------------
        else:
            print("  Regrow segments from watershed basin seeds")
            if regrow_by_propagating:
                shrunken = -1 * np.ones(len(depths))
                for iseed, seed_index in enumerate(seed_indices):
                    shrunken[seed_index] = iseed
            else:
                seed_lists = [[x] for x in seed_indices]

        #-----------------------------------------------------------------------
        # Growth by label propagation or segmentation
        #-----------------------------------------------------------------------
        if regrow_by_propagating:
            regrown = propagate(points, faces, folds, shrunken, folds,
                                max_iters=1000, tol=0.001, sigma=5)
        else:
            regrown = segment(indices_folds, neighbor_lists,
                              min_subfold_size, seed_lists)
        subfolds = np.array(subfolds)
        subfolds[regrown > -1] = regrown[regrown > -1] + np.max(folds) + 1

    #---------------------------------------------------------------------------
    # Renumber folds so they are sequential
    #---------------------------------------------------------------------------
    renumber_subfolds = -1 * np.ones(len(folds))
    subfold_numbers = [int(x) for x in np.unique(subfolds) if x > -1]
    for i_subfold, n_subfold in enumerate(subfold_numbers):
        subfold = [i for i,x in enumerate(subfolds) if x == n_subfold]
        renumber_subfolds[subfold] = i_subfold
    subfolds = renumber_subfolds
    n_subfolds = i_subfold + 1

    # Print statement
    print('  ...Extracted {0} subfolds ({1:.2f} seconds)'.
          format(n_subfolds, time() - t0))

    #---------------------------------------------------------------------------
    # Return subfolds, number of subfolds, file name
    #---------------------------------------------------------------------------
    if save_file:
        subfolds_file = os.path.join(os.getcwd(), 'subfolds.vtk')
        rewrite_scalars(depth_file, subfolds_file, subfolds, 'subfolds', subfolds)
    else:
        subfolds_file = None

    return subfolds.tolist(), n_subfolds, subfolds_file


#===============================================================================
# Normalize depths in each fold
#===============================================================================
def normalize_fold_depths(depth_file, folds_or_file, save_file=False):

    """
    Normalize depths in each fold.

    Parameters
    ----------
    depth_file : string
        name of VTK file with a depth value for each vertex
    folds_or_file : list or string
        fold number for each vertex or name of VTK file containing folds scalars
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    depth_folds : list of floats
        depth values for fold numbers >-1, normalized for each fold
    depth_folds_file : string (if save_file)
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
    >>> depth_folds_file = os.path.join(os.getcwd(),
    >>>     os.path.basename(depth_file).strip('vtk') + 'norm.vtk')
    >>> os.system('mayavi2 -m Surface -d ' + depth_folds_file + '&')

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

    # Load fold numbers if folds_or_file is a string
    if isinstance(folds_or_file, str):
        folds, name = read_scalars(folds_or_file)
    elif isinstance(folds_or_file, list):
        folds = folds_or_file

    depths, name = read_scalars(depth_file, True, True)

    depth_folds = -1 * np.ones(len(folds))

    unique_fold_IDs = np.unique(folds)
    unique_fold_IDs = [x for x in unique_fold_IDs if x >= 0]

    for fold_ID in unique_fold_IDs:
        indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
        if indices_fold:

            # Normalize fold depth values
            depth_folds[indices_fold] = depths[indices_fold] / np.max(depths[indices_fold])

    #---------------------------------------------------------------------------
    # Return folds with normalized depth, number of folds, file name
    #---------------------------------------------------------------------------
    if save_file:

        depth_folds_file = os.path.join(os.getcwd(),
            os.path.basename(depth_file).strip('vtk') + 'norm.folds.vtk')
        rewrite_scalars(depth_file, depth_folds_file,
                        depth_folds, 'depth_folds', depth_folds)

    else:
        depth_folds_file = None

    return depth_folds.tolist(), depth_folds_file
