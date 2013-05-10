#!/usr/bin/env python
"""
Miscellaneous computations.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def sigmoid(values, gain, shift):
    """
    Map values with sigmoid function to range [0,1].

    Y(t) = 1/(1 + exp(-gain*(values - shift))
    """
    import numpy as np

    tiny = 0.000000001

    # Make sure argument is a numpy array
    if type(values) != np.ndarray:
        values = np.array(values)

    return 1.0 / (1.0 + np.exp(-gain * (values - shift)) + tiny)

def find_segment_endpoints(indices, neighbor_lists, likelihoods, step=1):
    """
    Find endpoints in a region of connected vertices.

    These points are intended to serve as endpoints of fundus curves
    running along high-likelihood paths within a region (fold).
    This algorithm iteratively propagates from an initial high-likelihood seed
    toward the boundary of a region within thresholded subregions of
    decreasing likelihood. The first boundary point that is reached for each
    segmented branch serves as an endpoint.

    Note ::

        This algorithm suffers from the following problems:
        1.  If multiple boundary points are reached simultaneously,
            the choice of highest likelihood among them might not be appropriate.
        2.  The boundary may be reached before any other branches are
            discovered, especially if other branches are separated by
            low likelihood (shallow) vertices.
        3.  Segmentation may reach the top of a fold's wall before reaching
            the tip of its branch.

    Steps ::

        Initialize:
            R: Region/remaining vertices to segment (initially fold vertices)
            P: Previous segment (initially the maximum likelihood point)
            X: Excluded segment
            N: New/next segment (initially P)
            E: indices to endpoint vertices (initially empty)

        For each decreasing threshold, run recursive function creep():

            Propagate P into R, and call these new vertices N.
            Propagate X into P, R, and N.
            Optionally remove points from N and R that are also in the expanded X.
            Remove P and N from R.
            Reassign P to X.
            If N is empty:
                Choose highest likelihood point in P as endpoint.
                Return endpoints E and remaining vertices R.
            else:
                Identify N_i different segments of N.
                For each segment N_i:
                    If N_i large enough or if max(i)==1:
                        Call recursive function creep() with new arguments.
                Return endpoints E and R, P, X, and N.

    Parameters
    ----------
    indices : list of integers
        indices of the vertices to segment (such as a fold in a surface mesh)
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    likelihoods : numpy array of floats
        fundus likelihood values for all vertices
    step : integer
        number of segmentation steps before assessing segments

    Returns
    -------
    indices_endpoints : list of integers
        indices of surface mesh vertices that are endpoints

    Examples
    --------
    >>> # Setup:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.fundi import find_endpoints
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> likelihood_file = os.path.join(path, 'arno', 'features', 'likelihoods.vtk')
    >>> likelihoods, name = read_scalars(likelihood_file, True, True)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> neighbor_lists = find_neighbors_from_file(depth_file)
    >>> step = 1
    >>> min_size = 50
    >>> #
    >>> #-----------------------------------------------------------------------
    >>> # Find endpoints on a single fold:
    >>> #-----------------------------------------------------------------------
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != -1]
    >>> indices_endpoints = find_endpoints(indices, neighbor_lists,
    >>>                                    likelihoods, step)
    >>> # Write results to VTK file and view:
    >>> likelihoods[indices_endpoints] = max(likelihoods) + 0.1
    >>> rewrite_scalars(likelihood_file, 'find_endpoints.vtk',
    >>>                 likelihoods, 'endpoints_on_likelihoods_in_fold', fold)
    >>> plot_vtk('find_endpoints.vtk')
    >>> #
    >>> #-----------------------------------------------------------------------
    >>> # Find endpoints on every fold in a hemisphere:
    >>> #-----------------------------------------------------------------------
    >>> plot_each_fold = False
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> fold_numbers = [x for x in np.unique(folds) if x != -1]
    >>> nfolds = len(fold_numbers)
    >>> endpoints = []
    >>> for ifold, fold_number in enumerate(fold_numbers):
    >>>     print('Fold {0} ({1} of {2})'.format(int(fold_number), ifold+1, nfolds))
    >>>     indices = [i for i,x in enumerate(folds) if x == fold_number]
    >>>     if len(indices) > min_size:
    >>>         indices_endpoints = find_endpoints(indices, neighbor_lists, likelihoods, step)
    >>>         endpoints.extend(indices_endpoints)
    >>>         # Plot each fold:
    >>>         if plot_each_fold:
    >>>             fold = -1 * np.ones(len(likelihoods))
    >>>             fold[indices] = 1
    >>>             likelihoods[indices_endpoints] = max(likelihoods) + 0.1
    >>>             rewrite_scalars(likelihood_file, 'find_endpoints.vtk',
    >>>                     likelihoods, 'endpoints_on_likelihoods_in_fold', fold)
    >>>             plot_vtk('find_endpoints.vtk')
    >>> E = -1 * np.ones(len(likelihoods))
    >>> E[endpoints] = 1
    >>> #
    >>> # Write results to VTK file and view:
    >>> rewrite_scalars(folds_file, 'find_endpoints.vtk',
    >>>                 E, 'endpoints_on_folds', folds)
    >>> plot_vtk('find_endpoints.vtk')

    """
    import numpy as np

    from mindboggle.labels.labels import extract_borders

    # Make sure arguments are numpy arrays
    if isinstance(likelihoods, list):
        likelihoods = np.array(likelihoods)

    # Parameters:
    min_size = 1
    xstep = 1

    # Threshold parameters:
    use_thresholds = True
    threshold_factor = 0.9
    min_threshold = 0.1

    # Recursive function for segmenting and finding endpoints:
    def creep(R, P, X, E, L, B, step, neighbor_lists, min_size=1):
        """
        Recursively segment a mesh, creeping toward its edges to find endpoints.

        Steps ::

            Propagate P into R, and call these new vertices N.
            Propagate X into P, R, and N.
            Remove points from N and R that are also in the expanded X.
            Remove P and N from R.
            Reassign P to X.
            If N is empty:
                Choose highest likelihood point in P as endpoint.
                Return endpoints E and remaining vertices R.
            else:
                Identify N_i different segments of N.
                For each segment N_i:
                    If N_i large enough or if max(i)==1:
                        Call recursive function creep() with new arguments.
                Return E, R, P, X, and N.

        Parameters
        ----------
        R : list of integers
            indices of vertices to segment (such as a fold)
        P : list of integers
            indices of previous segment vertices
        X : list of integers
            indices of segmented vertices to exclude from endpoint selection
        E: list of integers
            indices to endpoint vertices
        L : numpy array of floats
            likelihood values for all vertices
        step : integer
            number of segmentation steps before assessing segments
        neighbor_lists : list of lists of integers
            indices to neighboring vertices for each vertex
        min_size : integer
            minimum number of vertices for an endpoint segment

        Returns
        -------
        R : list of integers
            remaining vertices to segment
        P : list of integers
            previous segment
        X : list of integers
            excluded segment
        E: list of integers
            endpoints

        """
        import numpy as np

        from mindboggle.utils.segment import segment

        # Expand X and exclude endpoint selection?:
        rmX = False

        #-----------------------------------------------------------------------
        # Propagate P into R, and call these new vertices N:
        #-----------------------------------------------------------------------
        PintoR = segment(R, neighbor_lists, min_region_size=1, seed_lists=[P],
                         keep_seeding=False, spread_within_labels=False,
                         labels=[], label_lists=[], values=[], max_steps=step)
        PN = [i for i,x in enumerate(PintoR) if x != -1]
        # Remove P (seeds) from N:
        N = list(frozenset(PN).difference(P))
        #print('  {0} vertices in the new segment'.format(len(N)))

        #-----------------------------------------------------------------------
        # Propagate X into R (including P and N):
        #-----------------------------------------------------------------------
        if rmX:
            if X:
                RPN = R[:]
                RPN.extend(PN)
                XintoR = segment(RPN, neighbor_lists, min_region_size=1,
                                 seed_lists=[X], keep_seeding=False,
                                 spread_within_labels=False, labels=[],
                                 label_lists=[], values=[], max_steps=xstep)
                X = [i for i,x in enumerate(XintoR) if x != -1]
                print('  {0} vertices spread from previously segmented'.format(len(X)))

                # Remove points from N and R that are also in the expanded X:
                N = list(frozenset(N).difference(X))
                R = list(frozenset(R).difference(X))

            # Reassign P to X:
            X.extend(P)

        # Remove P and N from R:
        R = list(frozenset(R).difference(P))
        R = list(frozenset(R).difference(N))

        #-----------------------------------------------------------------------
        # If N is empty, return endpoints:
        #-----------------------------------------------------------------------
        BandN = list(frozenset(B).intersection(N))

        if not N:
            pass

        elif BandN:

            # Choose highest likelihood point in P as endpoint:
            E.append(BandN[np.argmax(L[BandN])])

        #-----------------------------------------------------------------------
        # If N is not empty, assign as P and continue segmenting recursively:
        #-----------------------------------------------------------------------
        else:

            # Identify N_i different segments of N:
            N_segments = segment(N, neighbor_lists, min_region_size=1)
            unique_N = [x for x in np.unique(N_segments) if x!=-1]
            n_segments = len(unique_N)

            # For each segment N_i:
            for n in unique_N:
                N_i = [i for i,x in enumerate(N_segments) if x==n]

                # If N_i large enough or if max(i)==1:
                if len(N_i) >= min_size or n_segments==1:

                    # Call creep() with new arguments:
                    R, P, X, E = creep(R, N_i, X, E, L, B, step,
                                       neighbor_lists, min_size)

        # Return endpoints E and remaining vertices R:
        return R, P, X, E


    # Extract boundary:
    D = np.ones(len(likelihoods))
    D[indices] = 2
    B, foo1, foo2 = extract_borders(range(len(likelihoods)), D, neighbor_lists)

    # Initialize R, X, and E:
    R = []
    X = []
    E = []
    indices_endpoints = []

    # Initialize P and threshold with the maximum likelihood point:
    L = likelihoods
    Imax = indices[np.argmax(L[indices])]
    P = [Imax]
    threshold = L[Imax]

    # Include new vertices with lower likelihood values:
    if use_thresholds:

        # Iterate endpoint extraction until all vertices have been segmented:
        continue_loop = True
        while continue_loop:
            prev_threshold = threshold

            # If threshold above minimum, update R based on the threshold:
            if threshold > min_threshold:
                #if X:  threshold = threshold_factor * np.mean(L[X])
                threshold = threshold_factor * threshold
                T = [x for x in indices if L[x] >= threshold
                     if L[x] < prev_threshold]
                if not T:
                    decrease_threshold = True
                    while decrease_threshold:
                        threshold *= threshold_factor
                        T = [x for x in indices if L[x] >= threshold
                             if L[x] < prev_threshold]
                        if T or threshold < min_threshold:
                            decrease_threshold = False
                R.extend(T)

            # If threshold below minimum, update and exit:
            else:
                T = [x for x in indices if L[x] < prev_threshold]
                R.extend(T)
                continue_loop = False

            # Run recursive function creep() to return endpoints:
            R, P, X, E = creep(R, P, X, E, L, B, step, neighbor_lists, min_size)
            E = np.unique(E).tolist()

            # Print message:
            if len(R) == 1:
                str1 = 'vertex'
            else:
                str1 = 'vertices'
            if len(E) == 1:
                str2 = 'endpoint'
            else:
                str2 = 'endpoints'
            print('  {0} remaining {1}, {2} {3} (threshold: {4:0.3f})'.
                  format(len(R), str1, len(E), str2, threshold))

    # Don't use thresholds -- include all vertices:
    else:

        R = indices
        print('  Segment {0} vertices'.format(len(R)))

        # Run recursive function creep() to return endpoints:
        R, P, X, E = creep(R, P, X, E, L, B, step, neighbor_lists, min_size)

    indices_endpoints = E

    return indices_endpoints

def track(R, P, T, L, B, neighbor_lists):
    """
    Recursively run tracks along a mesh, through vertices of high likelihood.
    At each vertex, continue, branch, or terminate.

    Steps ::

        R is the set of remaining (untracked) vertices.
        Find the neighborhood N for point P in R.
        Remove N from R.
        For each neighborhood vertex N_i:
            Remove N_i from N.
            Find the neighbors for N_i also in N.
            If N_i has the maximum value in its restricted neighborhood:
                Call recursive function track() with N_i as P if N_i not in B.

    Parameters
    ----------
    R : list of integers
        indices of vertices (such as a fold in a surface mesh)
    P : integer
        index to vertex
    T : list of lists of pairs of integers
        index pairs are track segments
    L : numpy array of floats
        likelihood values for all vertices
    B : list of integers
        indices of boundary vertices for R
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex

    Returns
    -------
    R : list of integers
        remaining vertices
    T : list of lists of pairs of integers
        track segments

    """
    import numpy as np


    # Find the neighborhood N for point P in R:
    N = neighbor_lists[P]
    N = list(frozenset(N).intersection(R))
    print('N', N)
    if N:

        # Remove N from R:
        R = list(frozenset(R).difference(N))

        # For each neighborhood vertex N_i:
        Imax = np.argmax(L[N])
        print(Imax)
        for N_i in [N[Imax]]:

            # Find the neighbors of N_i also in N:
            N2 = list(frozenset(neighbor_lists[N_i]).intersection(N))
            print('N2', N2)
            if N2:

                # If N_i has the maximum value in its restricted neighborhood:
                if L[N_i] >= max(L[N2]):

                    # Add track segment:
                    T.append([P, N_i])
                    print('T', T)

            # Call recursive function track() with N_i as P if N_i not in B:
            if N_i not in B:
                R, T = track(R, N_i, T, L, B, neighbor_lists)

    return R, T

#-----------------------------------------------------------------------------
# Shrink segments
#-----------------------------------------------------------------------------
def shrink_segments(regions, segments, depths, shrink_factor=0.25,
                    only_multiple_segments=False):
    """
    Shrink segments in a segmented surface mesh by a fraction of its maximum
    depth, for all segments or for segments in regions with multiple segments.

    Parameters
    ----------
    regions : list or array of integers
        region IDs for all vertices, indicating inclusion in a region (default -1)
    segments : numpy array of integers
        segment IDs for all vertices, indicating inclusion in a segment (default -1)
    depths : numpy array of floats
        depth values for all vertices (default -1)
    shrink_factor : float
        shrink each region of connected vertices to this fraction
        of its maximum depth
    only_multiple_segments : Boolean
        shrink only segments in regions with multiple segments
        (otherwise shrink all segments)

    Returns
    -------
    shrunken_segments : list of integers
        shrunken segment numbers for all vertices (default -1)
        -- non-shrunken segments are removed

    Examples
    --------
    >>> # Segment folds with watershed(), then shrink these segments:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.utils.segment import watershed
    >>> from mindboggle.x.misc import shrink_segments
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> indices = np.where(depths > 0.11)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> segments = watershed(depths, points, indices, neighbor_lists, min_size=1,
    >>>     depth_factor=0.25, depth_ratio=0.1, tolerance=0.01)[0]
    >>> #
    >>> shrink_factor = 0.25
    >>> shrunken_segments = shrink_segments(folds, segments, depths,
    >>>     shrink_factor, only_multiple_segments=True)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'shrink_segments.vtk',
    >>>     shrunken_segments, 'shrunken_segments', shrunken_segments)
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('shrink_segments.vtk')

    """
    import numpy as np

    print('Shrink segments')

    remove_fraction = 1 - shrink_factor
    shrunken_segments = -1 * np.ones(len(depths))

    # Make sure arguments are numpy arrays
    if not isinstance(segments, np.ndarray):
        segments = np.array(segments)
    if not isinstance(depths, np.ndarray):
        depths = np.array(depths)

    # Shrink only segments in regions with multiple segments
    if only_multiple_segments:
        print('  Shrink each segment to {0:.2f} of its depth for regions with '
              'multiple segments'.format(shrink_factor))

        # For each region
        unique_regions = [x for x in np.unique(regions) if x > -1]
        for n_region in unique_regions:

            # Check to see if there are multiple segments in the region
            indices_region = [i for i,x in enumerate(regions) if x == n_region]
            segments_in_region = [x for x in np.unique(segments[indices_region])
                                  if x > -1]
            if len(segments_in_region) > 1:

                # Shrink each segment in the region
                for n_segment in segments_in_region:
                    indices_segment = [i for i,x in enumerate(segments)
                                       if x == n_segment]
                    indices_segment = list(frozenset(indices_segment).intersection(indices_region))
                    depth_threshold = remove_fraction * np.max(depths[indices_segment])
                    indices_segment = [x for x in indices_segment
                                       if depths[x] > depth_threshold]
                    shrunken_segments[indices_segment] = n_segment

    # Shrink all segments
    else:
        print('  Shrink each segment to {0:.2f} of its depth'.format(shrink_factor))
        unique_segments = [x for x in np.unique(segments) if x != -1]
        for n_segment in unique_segments:
            indices_segment = [i for i,x in enumerate(segments) if x == n_segment]
            depth_threshold = remove_fraction * np.max(depths[indices_segment])
            indices_segment = [x for x in indices_segment
                               if depths[x] > depth_threshold]
            shrunken_segments[indices_segment] = n_segment

    return shrunken_segments
