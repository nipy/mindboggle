#!/usr/bin/env python
"""
Functions to extract folds.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

<<<<<<< HEAD
#=============================================================================
# Extract folds
#=============================================================================
def extract_folds(depth_file, min_vertices=10000, min_fold_size=50, 
                  do_fill_holes=False, min_hole_depth=0.001, 
                  save_file=False):
=======
#===============================================================================
# Extract folds
#===============================================================================
def extract_folds(depth_file, neighbor_lists=[], min_fold_size=1,
                  extract_subfolds=True, save_file=False):
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    """
    Use depth to extract folds from a triangular surface mesh.

    Steps ::
<<<<<<< HEAD
        1. Compute histogram of depth measures.
        2. Define a depth threshold and find the deepest vertices.
        3. Segment deep vertices as an initial set of folds.
        4. Remove small folds.
        5. Find and fill holes in the folds (optional).
        6. Renumber folds.

    Step 2 ::
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
        To extract an initial set of deep vertices from the surface mesh,
        we anticipate that there will be a rapidly decreasing distribution
        of low depth values (on the outer surface) with a long tail
        of higher depth values (in the folds), so we smooth the histogram's
<<<<<<< HEAD
        bin values, convolve to compute slopes, and find the depth value
        for the first bin with slope = 0. This is our threshold.

    Step 5 ::
        The folds could have holes in areas shallower than the depth threshold.
        Calling fill_holes() could accidentally include very shallow areas
        (in an annulus-shaped fold, for example), so we include the argument
        exclude_range to check for any values from zero to min_hole_depth;
        holes are not filled if they contains values within this range.
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Parameters
    ----------
    depth_file : string
        surface mesh file in VTK format with faces and depth scalar values
<<<<<<< HEAD
    min_fold_size : integer
        minimum fold size (number of vertices)
    do_fill_holes : Boolean
        fill holes in the folds?
    min_hole_depth : float
        largest non-zero depth value that will stop a hole from being filled
=======
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
        -- if empty list, construct from depth_file
    min_fold_size : int
        minimum fold size (number of vertices)
    extract_subfolds : Boolean
        segment folds into subfolds?
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    save_file : Boolean
        save output VTK file?

    Returns
    -------
<<<<<<< HEAD
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
=======
    folds : array of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    n_folds :  int
        number of folds
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    folds_file : string (if save_file)
        name of output VTK file with fold IDs (-1 for non-fold vertices)

    Examples
    --------
    >>> import os
<<<<<<< HEAD
    >>> import numpy as np
    >>> import pylab
    >>> from scipy.ndimage.filters import gaussian_filter1d
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> from mindboggle.features.folds import extract_folds
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = 'travel_depth.vtk' #os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> min_vertices = 10000
    >>> min_fold_size = 50
    >>> do_fill_holes = False #True
    >>> min_hole_depth = 0.001
    >>> save_file = True
    >>> #
    >>> folds, n_folds, thr, bins, bin_edges, folds_file = extract_folds(depth_file,
    >>>     min_vertices, min_fold_size, do_fill_holes, min_hole_depth, save_file)
    >>> #
    >>> # View folds:
    >>> plot_surfaces('folds.vtk')
    >>> # Plot histogram and depth threshold:
    >>> depths, name = read_scalars(depth_file)
    >>> nbins = np.round(len(depths) / 100.0)
    >>> a,b,c = pylab.hist(depths, bins=nbins)
    >>> pylab.plot(thr*np.ones((100,1)), np.linspace(0, max(bins), 100), 'r.')
    >>> pylab.show()
    >>> # Plot smoothed histogram:
    >>> bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    >>> pylab.plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')
    >>> pylab.show()

    """
    import os
    import sys
=======
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
    import os
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    import numpy as np
    from time import time
    from scipy.ndimage.filters import gaussian_filter1d
    from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk
<<<<<<< HEAD
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.utils.morph import fill_holes
    from mindboggle.utils.segment import segment
=======
    from mindboggle.utils.mesh import fill_holes
    from mindboggle.labels.segment import segment, watershed, shrink_segments
    if not len(neighbor_lists):
        from mindboggle.utils.mesh import find_neighbors

    do_fill_holes = True
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    print("Extract folds in surface mesh")
    t0 = time()

<<<<<<< HEAD
    #-------------------------------------------------------------------------
    # Load depth values for all vertices
    #-------------------------------------------------------------------------
    faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
        return_first=True, return_array=True)

    #-------------------------------------------------------------------------
    # Find neighbors for each vertex
    #-------------------------------------------------------------------------
    neighbor_lists = find_neighbors(faces, npoints)

    #-------------------------------------------------------------------------
    # Compute histogram of depth measures
    #-------------------------------------------------------------------------
    if npoints > min_vertices:
        nbins = np.round(npoints / 100.0)
    else:
        sys.err("  Expecting at least {0} vertices to create depth histogram".
            format(min_vertices))
    bins, bin_edges = np.histogram(depths, bins=nbins)

    #-------------------------------------------------------------------------
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    # Anticipating that there will be a rapidly decreasing distribution
    # of low depth values (on the outer surface) with a long tail of higher
    # depth values (in the folds), smooth the bin values (Gaussian), convolve
    # to compute slopes, and find the depth for the first bin with slope = 0.
<<<<<<< HEAD
    #-------------------------------------------------------------------------
    bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    window = [-1, 0, 1]
    bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
    ibins0 = np.where(bin_slopes == 0)[0]
    if ibins0.shape:
=======
    #---------------------------------------------------------------------------
    bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    #>>> # Plot smoothed histogram:
    #>>> plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')
    window = [-1, 0, 1]
    bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
    ibins0 = np.where(bin_slopes == 0)[0]
    if ibins0.size:
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
        depth_threshold = bin_edges[ibins0[0]]
    else:
        depth_threshold = np.median(depths)

<<<<<<< HEAD
    #-------------------------------------------------------------------------
    # Find the deepest vertices
    #-------------------------------------------------------------------------
    indices_deep = [i for i,x in enumerate(depths) if x >= depth_threshold]
    if indices_deep:

        #---------------------------------------------------------------------
        # Segment deep vertices as an initial set of folds
        #---------------------------------------------------------------------
=======
    #---------------------------------------------------------------------------
    # Find the deepest vertices
    #---------------------------------------------------------------------------
    indices_deep = [i for i,x in enumerate(depths) if x >= depth_threshold]
    if indices_deep:

        #-----------------------------------------------------------------------
        # Segment deep vertices as an initial set of folds
        #-----------------------------------------------------------------------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
        print("  Segment vertices deeper than {0:.2f} as folds".format(depth_threshold))
        t1 = time()
        folds = segment(indices_deep, neighbor_lists)
        # Slightly slower alternative -- fill boundaries:
        #regions = -1 * np.ones(len(points))
        #regions[indices_deep] = 1
<<<<<<< HEAD
        #folds = segment_by_filling_borders(regions, neighbor_lists)
        print('  ...Segmented folds ({0:.2f} seconds)'.format(time() - t1))

        #---------------------------------------------------------------------
        # Remove small folds
        #---------------------------------------------------------------------
        if min_fold_size > 1:
            print('  Remove folds smaller than {0}'.format(min_fold_size))
            unique_folds = [x for x in np.unique(folds) if x != -1]
            for nfold in unique_folds:
                indices_fold = [i for i,x in enumerate(folds) if x == nfold]
                if len(indices_fold) < min_fold_size:
                    folds[indices_fold] = -1

        #---------------------------------------------------------------------
        # Find and fill holes in the folds
        # Note: Surfaces surrounded by folds can be mistaken for holes,
        #       so exclude_range includes outer surface values close to zero.
        #---------------------------------------------------------------------
        if do_fill_holes:
            print("  Find and fill holes in the folds")
            folds = fill_holes(folds, neighbor_lists, values=depths,
                               exclude_range=[0, min_hole_depth])

        #---------------------------------------------------------------------
        # Renumber folds so they are sequential
        #---------------------------------------------------------------------
        renumber_folds = -1 * np.ones(len(folds))
        fold_numbers = [int(x) for x in np.unique(folds) if x != -1]
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
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

<<<<<<< HEAD
    folds = [int(x) for x in folds]

    #-------------------------------------------------------------------------
    # Return folds, number of folds, file name
    #-------------------------------------------------------------------------
=======
    #---------------------------------------------------------------------------
    # Return folds, number of folds, file name
    #---------------------------------------------------------------------------
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    if save_file:

        folds_file = os.path.join(os.getcwd(), 'folds.vtk')
        rewrite_scalars(depth_file, folds_file, folds, 'folds', folds)

<<<<<<< HEAD
        if not os.path.exists(folds_file):
            raise(IOError(folds_file + " not found"))

    else:
        folds_file = None

    return folds, n_folds, depth_threshold, bins, bin_edges, folds_file


#=============================================================================
# Extract subfolds
#=============================================================================
def extract_subfolds(depth_file, folds, min_size=10, depth_factor=0.25,
                     depth_ratio=0.1, tolerance=0.01, save_file=False):
    """
    Use depth to segment folds into subfolds in a triangular surface mesh.

    Note ::

        The function extract_sulci() performs about the same whether folds
        or subfolds are used as input.  The latter leads to some loss of
        small subfolds and possibly holes for small subfolds in the middle
        of other subfolds.

    Note about the watershed() function:
    The watershed() function performs individual seed growing from deep seeds,
    repeats segmentation from the resulting seeds until each seed's segment
    touches a boundary. The function segment() fills in the rest. Finally
    segments are joined if their seeds are too close to each other.
    Despite these precautions, the order of seed selection in segment() could
    possibly influence the resulting borders between adjoining segments.
    [The propagate() function is slower and insensitive to depth,
     but is not biased by seed order.]
=======
    else:
        folds_file = None

    return folds, n_folds, folds_file


#===============================================================================
# Normalize depths in each fold
#===============================================================================
def normalize_fold_depths(depth_file, folds, save_file=False):

    """
    Normalize depths in each fold.
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Parameters
    ----------
    depth_file : string
<<<<<<< HEAD
        surface mesh file in VTK format with faces and depth scalar values
    folds : list of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    min_size : integer
        minimum number of vertices for a subfold
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
=======
        name of VTK file with a depth value for each vertex
    folds : list
        fold ID for each vertex
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
    save_file : Boolean
        save output VTK file?

    Returns
    -------
<<<<<<< HEAD
    subfolds : list of integers
        fold numbers for all vertices (-1 for non-fold vertices)
    n_subfolds :  int
        number of subfolds
    subfolds_file : string (if save_file)
        name of output VTK file with fold IDs (-1 for non-fold vertices)
=======
    depth_folds : list of floats
        depth values for fold numbers >-1, normalized for each fold
    depth_folds_file : string (if save_file)
        name of output VTK file with normalized depth values in each fold
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    Examples
    --------
    >>> import os
<<<<<<< HEAD
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.folds import extract_subfolds
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> min_size = 10
    >>> depth_factor = 0.5
    >>> depth_ratio = 0.1
    >>> tolerance = 0.01
    >>> #
    >>> subfolds, n_subfolds, subfolds_file = extract_subfolds(depth_file,
    >>>     folds, min_size, depth_factor, depth_ratio, tolerance, True)
    >>> #
    >>> # View:
    >>> rewrite_scalars(depth_file, 'subfolds.vtk', subfolds, 'subfolds', subfolds)
    >>> plot_surfaces('subfolds.vtk')
=======
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
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

    """
    import os
    import numpy as np
<<<<<<< HEAD
    from time import time
    from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.utils.segment import segment, propagate, watershed

    print("Segment folds into subfolds")
    t0 = time()

    #-------------------------------------------------------------------------
    # Load depth values for all vertices
    #-------------------------------------------------------------------------
    faces, lines, indices, points, npoints, depths, \
        name, input_vtk = read_vtk(depth_file, return_first=True, return_array=True)

    #-------------------------------------------------------------------------
    # Find neighbors for each vertex
    #-------------------------------------------------------------------------
    neighbor_lists = find_neighbors(faces, npoints)

    #-------------------------------------------------------------------------
    # Segment folds into "watershed basins"
    #-------------------------------------------------------------------------
    indices_folds = [i for i,x in enumerate(folds) if x != -1]
    subfolds, seed_indices = watershed(depths, points, indices_folds,
                                 neighbor_lists, min_size, depth_factor=0.25,
                                 depth_ratio=0.1, tolerance=0.01, regrow=True)

    # Print statement
    n_subfolds = len([x for x in np.unique(subfolds) if x != -1])
    print('  ...Extracted {0} subfolds ({1:.2f} seconds)'.
          format(n_subfolds, time() - t0))

    #-------------------------------------------------------------------------
    # Return subfolds, number of subfolds, file name
    #-------------------------------------------------------------------------
    if save_file:
        subfolds_file = os.path.join(os.getcwd(), 'subfolds.vtk')
        rewrite_scalars(depth_file, subfolds_file, subfolds, 'subfolds', subfolds)

        if not os.path.exists(subfolds_file):
            raise(IOError(subfolds_file + " not found"))

    else:
        subfolds_file = None

    return subfolds, n_subfolds, subfolds_file
=======
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

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

    return depth_folds, depth_folds_file
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4
