#!/usr/bin/env python
"""
Functions for segmenting a surface mesh.

Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def propagate(points, faces, region, seeds, labels,
              max_iters=500, tol=0.001, sigma=10, background_value=-1,
              verbose=False):
    """
    Propagate labels to segment surface into contiguous regions,
    starting from seed vertices.

    Imports : mindboggle.guts.rebound

    Parameters
    ----------
    points : array (or list) of lists of three integers
        coordinates for all vertices
    faces : list of lists of three integers
        indices to three vertices per face (indices start from zero)
    region : list (or array) of integers
        values > background_value: inclusion in a region for all vertices
    seeds : numpy array of integers
        seed numbers for all vertices
    labels : numpy array of integers
        label numbers for all vertices
    max_iters : integer
        maximum number of iterations to run graph-based learning algorithm
    tol: float
        threshold to assess convergence of the algorithm
    sigma: float
        gaussian kernel parameter
    background_value : integer
        background value
    verbose : bool
        print statements?

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Propagate labels between label boundary segments in a single fold:
    >>> import numpy as np
    >>> import mindboggle.guts.rebound as rb
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.guts.segment import extract_borders
    >>> from mindboggle.guts.segment import propagate
    >>> from mindboggle.mio.vtks import read_scalars, read_vtk
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> dkt = DKTprotocol()
    >>> folds, name = read_scalars(folds_file, True, True)
    >>> points, f1,f2, faces, labels, f3, npoints, f4 = read_vtk(label_file,
    ...     True, True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> background_value = -1
    >>> # Limit number of folds to speed up the test:
    >>> limit_folds = True
    >>> if limit_folds:
    ...     fold_numbers = [4, 7] #[4, 6]
    ...     indices_fold = [i for i,x in enumerate(folds) if x in fold_numbers]
    ...     i0 = [i for i,x in enumerate(folds) if x not in fold_numbers]
    ...     folds[i0] = background_value
    ... else:
    ...     indices_fold = range(len(values))
    >>> # Extract the boundary for this fold:
    >>> indices_borders, label_pairs, foo = extract_borders(indices_fold,
    ...     labels, neighbor_lists, [], True)
    >>> # Select boundary segments in the sulcus labeling protocol:
    >>> seeds = background_value * np.ones(npoints)
    >>> for ilist,label_pair_list in enumerate(dkt.sulcus_label_pair_lists):
    ...     I = [x for i,x in enumerate(indices_borders)
    ...          if np.sort(label_pairs[i]).tolist() in label_pair_list]
    ...     seeds[I] = ilist
    >>> verbose = False
    >>> region = folds
    >>> max_iters = 500
    >>> tol = 0.001
    >>> sigma = 10
    >>> segments = propagate(points, faces, region, seeds, labels,
    ...                      max_iters, tol, sigma, background_value, verbose)
    >>> np.unique(segments)[0:10]
    array([ -1.,   3.,  12.,  23.])
    >>> len_segments = [len(np.where(segments == x)[0])
    ...                 for x in np.unique(segments) if x != background_value]
    >>> len_segments[0:10]
    [1152, 388, 116]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> rewrite_scalars(label_file, 'propagate.vtk',
    ...                 segments, 'segments', folds, background_value) # doctest: +SKIP
    >>> plot_surfaces('propagate.vtk') # doctest: +SKIP

    """
    import numpy as np
    from mindboggle.guts.mesh import keep_faces
    import mindboggle.guts.kernels as kernels
    import mindboggle.guts.rebound as rebound

    # Make sure arguments are numpy arrays:
    if not isinstance(seeds, np.ndarray):
        seeds = np.array(seeds)
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)
    if not isinstance(points, np.ndarray):
        points = np.array(points)

    if points.size and faces:
        segments = background_value * np.ones(len(points))
        indices_region = [i for i,x in enumerate(region)
                          if x != background_value]
        if indices_region:
            local_indices_region = background_value * np.ones(labels.shape)
            local_indices_region[indices_region] = list(range(len(indices_region)))

            if verbose:
                n_sets = len(np.unique([x for x in seeds
                                        if x != background_value]))
                if n_sets == 1:
                    print('Segment {0} vertices from 1 set of seed vertices'.
                          format(len(indices_region)))
                else:
                    print('Segment {0} vertices from {1} sets of seed '
                          'vertices'.format(len(indices_region), n_sets))

            # Remove faces whose 3 vertices are not among specified indices:
            refaces = keep_faces(faces, indices_region)

            # Set up rebound Bounds class instance:
            B = rebound.Bounds()
            if refaces:
                B.Faces = np.array(refaces)
                B.Indices = local_indices_region
                B.Points = points[indices_region]
                B.Labels = labels[indices_region]
                B.seed_labels = seeds[indices_region]
                B.num_points = len(B.Points)
                B.verbose = verbose

                # Propagate seed IDs from seeds:
                B.graph_based_learning(method='propagate_labels',
                                       realign=False,
                                       kernel=kernels.rbf_kernel,
                                       sigma=sigma,
                                       max_iters=max_iters,
                                       tol=tol,
                                       vis=False,
                                       verbose=verbose)

                # Assign maximum probability seed IDs to each point of region:
                max_prob_labels = B.assign_max_prob_label(verbose=False)

                # Return segment IDs in original vertex array:
                segments[indices_region] = max_prob_labels
            else:
                if verbose:
                    print("  No faces")
    else:
        segments = []

    return segments


def segment_regions(vertices_to_segment, neighbor_lists, min_region_size=1,
                    seed_lists=[], keep_seeding=False,
                    spread_within_labels=False, labels=[], label_lists=[],
                    values=[], max_steps='', background_value=-1, verbose=False):
    """
    Segment vertices of surface into contiguous regions by seed growing,
    starting from zero or more lists of seed vertices.

    Parameters
    ----------
    vertices_to_segment : list of integers
        indices to mesh vertices to be segmented
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    min_region_size : integer
        minimum size of segmented set of vertices
    seed_lists : list of lists, or empty list
        each list contains indices to seed vertices to segment vertices_to_segment
    keep_seeding : bool
        grow from new seeds even after all seed lists have fully grown
    spread_within_labels : bool
        grow seeds only by vertices with labels in the seed labels?
    labels : list of integers (required only if spread_within_labels)
        label numbers for all vertices
    label_lists : list of lists of integers (required only if spread_within_labels)
        List of unique labels for each seed list to grow into
        (If empty, set to unique labels for each seed list)
    values : list of floats (default empty)
        values for all vertices for use in preferentially directed segmentation
        (segment in direction of lower values)
    max_steps : integer (or empty string for infinity)
        maximum number of segmentation steps to take for each seed list
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Segment deep regions with or without seeds:
    >>> import numpy as np
    >>> from mindboggle.guts.segment import segment_regions
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> background_value = -1
    >>> urls, fetch_data = prep_tests()
    >>> depth_file = fetch_data(urls['left_travel_depth'])
    >>> f1,f2,f3, faces, depths, f4, npoints, t5 = read_vtk(depth_file,
    ...                                                     True, True)
    >>> vertices_to_segment = np.where(depths > 0.50)[0].tolist()  # (sped up)
    >>> neighbor_lists = find_neighbors(faces, npoints)

    Example 1: without seed lists

    >>> segments = segment_regions(vertices_to_segment, neighbor_lists)
    >>> len(np.unique(segments))
    92
    >>> len_segments = [len(np.where(segments == x)[0])
    ...                 for x in np.unique(segments) if x != background_value]
    >>> len_segments[0:10]
    [26631, 110928, 4, 1399, 1274, 5, 139, 255, 12, 5]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> rewrite_scalars(depth_file, 'segment.vtk', segments, 'segments', [], -1) # doctest: +SKIP
    >>> plot_surfaces('segment.vtk') # doctest: +SKIP

    Example 2: with seed lists

    >>> from mindboggle.guts.segment import extract_borders
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> dkt = DKTprotocol()
    >>> label_lists = [np.unique(np.ravel(x)) for x in dkt.sulcus_label_pair_lists]
    >>> labels_file = fetch_data(urls['left_freesurfer_labels'])
    >>> labels, name = read_scalars(labels_file)
    >>> indices_borders, label_pairs, foo = extract_borders(vertices_to_segment,
    ...     labels, neighbor_lists, ignore_values=[], return_label_pairs=True)
    >>> seed_lists = []
    >>> for label_pair_list in dkt.sulcus_label_pair_lists:
    ...     seed_lists.append([x for i,x in enumerate(indices_borders)
    ...                        if np.sort(label_pairs[i]).tolist()
    ...                        in label_pair_list])
    >>> keep_seeding = True
    >>> spread_within_labels = True
    >>> values = []
    >>> max_steps = ''
    >>> background_value = -1
    >>> verbose = False
    >>> segments = segment_regions(vertices_to_segment,
    ...     neighbor_lists, 1, seed_lists, keep_seeding, spread_within_labels,
    ...     labels, label_lists, values, max_steps, background_value, verbose)
    >>> len(np.unique(segments))
    122
    >>> np.unique(segments)[0:10]
    array([ -1.,   1.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  11.])
    >>> len_segments = [len(np.where(segments == x)[0])
    ...                 for x in np.unique(segments) if x != background_value]
    >>> len_segments[0:10]
    [6962, 8033, 5965, 5598, 7412, 3636, 3070, 5244, 3972, 6144]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> rewrite_scalars(depth_file, 'segment_regions.vtk', segments,
    ...     'segments', [], -1) # doctest: +SKIP
    >>> plot_surfaces('segment_regions.vtk') # doctest: +SKIP

    """
    import numpy as np

    verbose = False

    # Make sure arguments are lists:
    if isinstance(vertices_to_segment, np.ndarray):
        vertices_to_segment = vertices_to_segment.tolist()
    if isinstance(labels, np.ndarray):
        labels = [int(x) for x in labels]
    if isinstance(values, np.ndarray):
        values = values.tolist()

    #-------------------------------------------------------------------------
    # If seed_lists is empty, select first vertex from vertices_to_segment
    # (single vertex selection does not affect result -- see below*):
    #-------------------------------------------------------------------------
    if seed_lists:
        select_single_seed = False
        if verbose:
            if len(seed_lists) == 1:
                print('    Segment {0} vertices from seed vertices'.
                      format(len(vertices_to_segment)))
            else:
                print('    Segment {0} vertices from {1} sets of seed vertices'.
                      format(len(vertices_to_segment), len(seed_lists)))
    else:
        select_single_seed = True
        seed_lists = [[vertices_to_segment[0]]]
        if verbose:
            print('    Segment {0} vertices from first vertex as initial seed'.
                  format(len(vertices_to_segment)))

    #-------------------------------------------------------------------------
    # Initialize variables, including the list of vertex indices for each region,
    # vertex indices for all regions, and Boolean list indicating which regions
    # are fully grown, number of segments, etc.:
    #-------------------------------------------------------------------------
    segments = background_value * np.ones(len(neighbor_lists))
    region_lists = [[] for x in seed_lists]
    all_regions = []
    fully_grown = [False for x in seed_lists]
    new_segment_index = 0
    counter = 0
    if isinstance(max_steps, str):
        max_steps = np.Inf

    #-------------------------------------------------------------------------
    # If label_lists empty, set to unique labels for each seed list:
    #-------------------------------------------------------------------------
    if spread_within_labels:
        if not len(label_lists):
            label_lists = []
            for seed_list in seed_lists:
                seed_labels = np.unique([labels[x] for x in seed_list])
                label_lists.append(seed_labels)

    #-------------------------------------------------------------------------
    # Loop until all of the seed lists have grown to their full extent:
    #-------------------------------------------------------------------------
    count = 0
    while not all(fully_grown):
        # Loop through seed lists over and over again:
        for ilist, seed_list in enumerate(seed_lists):
            # If seed list empty
            if not len(seed_list):
                fully_grown[ilist] = True
            # If seed list not fully grown:
            if not fully_grown[ilist]:

                # Add seeds to region:
                region_lists[ilist].extend(seed_list)
                all_regions.extend(seed_list)

                # Remove seeds from vertices to segment:
                vertices_to_segment = list(frozenset(vertices_to_segment).
                difference(seed_list))

                if vertices_to_segment:

                    # Find neighbors of each seed with lower values than the seed:
                    if values:
                        neighbors = []
                        for seed in seed_list:
                            seed_neighbors = neighbor_lists[seed]
                            seed_neighbors = [x for x in seed_neighbors
                                              if values[x] <= values[seed]]
                            if seed_neighbors:
                                neighbors.extend(seed_neighbors)
                    else:
                        # Find neighbors of seeds:
                        neighbors = []
                        [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                    # Select neighbors that have not been previously selected
                    # and are among the vertices to segment:
                    seed_list = list(frozenset(neighbors).intersection(vertices_to_segment))
                    seed_list = list(frozenset(seed_list).difference(all_regions))

                else:
                    seed_list = []

                # If there are seeds remaining:
                if seed_list and count < max_steps:

                    # Select neighbors with the same labels
                    # as the initial seed labels:
                    if spread_within_labels:
                        seed_list = [x for x in seed_list
                                     if labels[x] in label_lists[ilist]]

                    # Continue growing seed list:
                    seed_lists[ilist] = seed_list
                    count += 1

                # If there are no seeds remaining:
                else:
                    # Stop growing seed list (see exception below):
                    fully_grown[ilist] = True

                    # If the region size is large enough:
                    size_region = len(region_lists[ilist])
                    if size_region >= min_region_size:

                        # Assign ID to segmented region and increment ID:
                        if select_single_seed:
                            new_segment_index = counter
                            counter += 1
                        else:
                            new_segment_index = ilist
                        segments[region_lists[ilist]] = new_segment_index

                        # Display current number and size of region:
                        if verbose and size_region > 1:
                            if len(seed_lists) == 1 and vertices_to_segment:
                                print("      {0} vertices remain".
                                      format(len(vertices_to_segment)))
                            else:
                                print("      Region {0}: {1} vertices ({2} remain)".
                                      format(int(new_segment_index), size_region,
                                             len(vertices_to_segment)))

                    # If selecting a single seed, continue growing
                    # if there are more vertices to segment:
                    if select_single_seed and count < max_steps:
                        if len(vertices_to_segment) >= min_region_size:
                            fully_grown[0] = False
                            seed_lists[0] = [vertices_to_segment[0]]
                            region_lists[0] = []

    #-------------------------------------------------------------------------
    # Keep growing from new seeds even after all seed lists have fully grown:
    #-------------------------------------------------------------------------
    if keep_seeding and len(vertices_to_segment) >= min_region_size:
        if verbose:
            print('    Keep seeding to segment {0} remaining vertices'.
                  format(len(vertices_to_segment)))

        # Select first unsegmented vertex as new seed:
        seed_list = [vertices_to_segment[0]]

        # Loop until the seed list has grown to its full extent:
        new_segment_index = ilist + 1
        region = []
        while len(vertices_to_segment) >= min_region_size:

            # Add seeds to region:
            region.extend(seed_list)
            all_regions.extend(seed_list)

            # Remove seeds from vertices to segment:
            vertices_to_segment = list(frozenset(vertices_to_segment).
            difference(seed_list))
            if vertices_to_segment:

                # Identify neighbors of seeds:
                neighbors = []
                [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                # Select neighbors that have not been previously selected
                # and are among the vertices to segment:
                seed_list = list(frozenset(vertices_to_segment).intersection(neighbors))
                seed_list = list(frozenset(seed_list).difference(all_regions))
            else:
                seed_list = []

            # If there are no seeds remaining:
            if not len(seed_list):

                # If the region size is large enough:
                size_region = len(region)
                if size_region >= min_region_size:

                    # Assign ID to segmented region and increment ID:
                    segments[region] = new_segment_index
                    new_segment_index += 1

                    # Display current number and size of region:
                    if verbose and size_region > 1:
                        print("      {0} vertices remain".
                              format(len(vertices_to_segment)))

                # Select first unsegmented vertex as new seed:
                if len(vertices_to_segment) >= min_region_size:
                    seed_list = [vertices_to_segment[0]]
                    region = []

    return segments


def segment_by_filling_borders(regions, neighbor_lists, background_value=-1,
                               verbose=False):
    """
    Fill borders (contours) on a surface mesh
    to segment vertices into contiguous regions.

    Steps ::
        1. Extract region borders (assumed to be closed contours)
        2. Segment borders into separate, contiguous borders
        3. For each boundary
            4. Find the neighbors to either side of the boundary
            5. Segment the neighbors into exterior and interior sets of neighbors
            6. Find the interior (smaller) sets of neighbors
            7. Fill the contours formed by the interior neighbors

    Parameters
    ----------
    regions : numpy array of integers
        region numbers for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Segment folds by extracting their borders and filling them in separately:
    >>> import numpy as np
    >>> from mindboggle.guts.segment import segment_by_filling_borders
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> depth_file = fetch_data(urls['left_travel_depth'])
    >>> f1,f2,f3, faces, depths, f4, npoints, input_vtk = read_vtk(depth_file,
    ...                                                            True, True)
    >>> background_value = -1
    >>> regions = background_value * np.ones(npoints)
    >>> regions[depths > 0.50] = 1
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> verbose = False
    >>> segments = segment_by_filling_borders(regions, neighbor_lists,
    ...                                       background_value, verbose)
    >>> len(np.unique(segments))
    19
    >>> len_segments = []
    >>> for useg in np.unique(segments):
    ...     len_segments.append(len(np.where(segments == useg)[0]))
    >>> len_segments[0:10]
    [19446, 8619, 13846, 23, 244, 101687, 16, 792, 210, 76]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> rewrite_scalars(depth_file, 'segment_by_filling_borders.vtk',
    ...                 segments, 'segments', [], -1) # doctest: +SKIP
    >>> plot_surfaces('segment_by_filling_borders.vtk') # doctest: +SKIP

    """
    import numpy as np
    from mindboggle.guts.segment import extract_borders
    from mindboggle.guts.segment import segment_regions

    include_boundary = False

    # Make sure arguments are numpy arrays
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    if verbose:
        print('Segment vertices using region borders')

    # Extract region borders (assumed to be closed contours)
    if verbose:
        print('  Extract region borders (assumed to be closed contours)')
    indices_borders, foo1, foo2 = extract_borders(list(range(len(regions))),
                                                  regions, neighbor_lists)
    # Extract background
    indices_background = list(frozenset(list(range(len(regions)))).
    difference(indices_borders))

    # Segment borders into separate, contiguous borders
    if verbose:
        print('  Segment borders into separate, contiguous borders')
    borders = segment_regions(indices_borders, neighbor_lists, 1, [], False,
                              False, [], [], [], '', background_value,
                              verbose)

    # For each boundary
    unique_borders = [x for x in np.unique(borders) if x != background_value]
    segments = background_value * np.ones(len(regions))
    for boundary_number in unique_borders:

        if verbose:
            print('  Boundary {0} of {1}:'.format(int(boundary_number),
                                                  len(unique_borders)))
        border_indices = [i for i,x in enumerate(borders)
                          if x == boundary_number]
        # Find the neighbors to either side of the boundary
        indices_neighbors = []
        [indices_neighbors.extend(neighbor_lists[i]) for i in border_indices]
        #indices_neighbors2 = indices_neighbors[:]
        #[indices_neighbors2.extend(neighbor_lists[i]) for i in indices_neighbors]
        indices_neighbors = list(frozenset(indices_neighbors).
        difference(indices_borders))

        # Segment the neighbors into exterior and interior sets of neighbors
        if verbose:
            print('    Segment the neighbors into exterior and interior sets '
                  'of neighbors')
        neighbors = segment_regions(indices_neighbors, neighbor_lists, 1, [],
                                    False, False, [], [], [], '',
                                    background_value, verbose)

        # Find the interior (smaller) sets of neighbors
        if verbose:
            print('    Find the interior (smaller) sets of neighbors')
        seed_lists = []
        unique_neighbors = [x for x in np.unique(neighbors)
                            if x != background_value]
        max_neighbor = 0
        max_len = 0
        for ineighbor, neighbor in enumerate(unique_neighbors):
            indices_neighbor = [i for i,x in enumerate(neighbors)
                                if x == neighbor]
            seed_lists.append(indices_neighbor)
            if len(indices_neighbor) > max_len:
                max_len = len(indices_neighbor)
                max_neighbor = ineighbor
        seed_lists = [x for i,x in enumerate(seed_lists) if i != max_neighbor]
        seed_list = []
        [seed_list.extend(x) for x in seed_lists if len(x) > 2]

        # Fill the contours formed by the interior neighbors
        if verbose:
            print('    Fill the contour formed by the interior neighbors')
        vertices_to_segment = list(frozenset(indices_background).
        difference(indices_borders))
        segmented_regions = segment_regions(vertices_to_segment,
                                            neighbor_lists, 1, [seed_list],
                                            False, False, [], [], [], '',
                                            background_value, verbose)

        if include_boundary:
            segmented_regions[border_indices] = 1

        segments[segmented_regions != background_value] = boundary_number

    return segments


def segment_rings(region, seeds, neighbor_lists, step=1, background_value=-1,
                  verbose=False):
    """
    Iteratively segment a region of surface mesh as concentric segments.

    Parameters
    ----------
    region : list of integers
        indices of region vertices to segment (such as a fold)
    seeds : list of integers
        indices of seed vertices
    neighbor_lists : list of lists of integers
        indices to neighboring vertices for each vertex
    step : integer
        number of segmentation steps before assessing segments
    background_value : integer
        background value
    verbose : bool
        print statements?

    Returns
    -------
    segments : list of lists of integers
        indices to vertices for each concentric segment

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.guts.mesh import find_neighbors_from_file
    >>> from mindboggle.guts.segment import extract_borders
    >>> from mindboggle.guts.segment import segment_rings
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_travel_depth'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> values, name = read_scalars(vtk_file, True, True)
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> background_value = -1
    >>> fold, name = read_scalars(folds_file)
    >>> indices = [i for i,x in enumerate(fold) if x != background_value]
    >>> # Initialize seeds with the boundary of thresholded indices:
    >>> use_threshold = True
    >>> if use_threshold:
    ...     # Threshold at the median depth or within maximum values in boundary:
    ...     threshold = np.median(values[indices]) #+ np.std(values[indices])
    ...     indices_high = [x for x in indices if values[x] >= threshold]
    ...     # Make sure threshold is within the maximum values of the boundary:
    ...     B = np.ones(len(values))
    ...     B[indices] = 2
    ...     borders, foo1, foo2 = extract_borders(list(range(len(B))), B, neighbor_lists)
    ...     borders = [x for x in borders if values[x] != background_value]
    ...     if list(frozenset(indices_high).intersection(borders)):
    ...         threshold = np.max(values[borders]) + np.std(values[borders])
    ...         indices_high = [x for x in indices if values[x] >= threshold]
    ...     # Extract threshold boundary vertices as seeds:
    ...     B = background_value * np.ones(len(values))
    ...     B[indices_high] = 2
    ...     seeds, foo1, foo2 = extract_borders(list(range(len(values))), B, neighbor_lists)
    ... # Or initialize P with the maximum value point:
    ... else:
    ...     seeds = [indices[np.argmax(values[indices])]]
    ...     indices_high = []
    >>> indices = list(frozenset(indices).difference(indices_high))
    >>> indices = list(frozenset(indices).difference(seeds))
    >>> step = 1
    >>> verbose = False
    >>> segments = segment_rings(indices, seeds, neighbor_lists, step,
    ...                          background_value, verbose)
    >>> len(segments)
    56
    >>> [len(x) for x in segments][0:10]
    [5540, 5849, 6138, 5997, 4883, 3021, 1809, 1165, 842, 661]
    >>> segments[0][0:10]
    [65539, 65540, 98308, 98316, 131112, 131121, 131122, 131171, 131175, 131185]

    Write results to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import read_scalars, rewrite_scalars # doctest: +SKIP
    >>> S = background_value * np.ones(len(values)) # doctest: +SKIP
    >>> for i, segment in enumerate(segments): S[segment] = i # doctest: +SKIP
    >>> rewrite_scalars(vtk_file, 'segment_rings.vtk', S, 'segment_rings',
    ...                 [], -1) # doctest: +SKIP
    >>> plot_surfaces('segment_rings.vtk') # doctest: +SKIP

    """
    from mindboggle.guts.segment import segment_regions

    segments = []
    while seeds:

        # Segment step-wise starting from seeds and through the region:
        seeds_plus_new = segment_regions(region, neighbor_lists, 1, [seeds],
                                         False, False, [], [], [],
                                         step, background_value, verbose)

        seeds_plus_new = [i for i,x in enumerate(seeds_plus_new)
                          if x != background_value]

        # Store the new segment after removing the previous segment:
        region = list(frozenset(region).difference(seeds))
        seeds = list(frozenset(seeds_plus_new).difference(seeds))
        if seeds:

            # Add the new segment and remove it from the region:
            segments.append(seeds)
            region = list(frozenset(region).difference(seeds))

    return segments

#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
def watershed(depths, points, indices, neighbor_lists, min_size=1,
              depth_factor=0.25, depth_ratio=0.1, tolerance=0.01, regrow=True,
              background_value=-1, verbose=False):
    """
    Segment surface vertices into contiguous regions by a watershed algorithm.

    Segment surface mesh into contiguous "watershed basins"
    by seed growing from an iterative selection of the deepest vertices.

    Steps ::

        1. Grow segments from an iterative selection of the deepest seeds.
        2. Regrow segments from the resulting seeds, until each seed's
            segment touches a boundary.
        3. Use the segment() function to fill in the rest.
        4. Merge segments if their seeds are too close to each other
            or their depths are very different.

    Note ::

        Despite the above precautions, the order of seed selection in segment()
        could possibly influence the resulting borders between adjoining
        segments (vs. propagate(), which is slower and insensitive to depth,
        but is not biased by seed order).

        In the example below, white spots indicate incomplete segmentation.

    Parameters
    ----------
    depths : numpy array of floats
        depth values for all vertices
    points : list of lists of floats
        each element is a list of 3-D coordinates of a vertex on a surface mesh
    indices : list of integers
        indices to mesh vertices to be segmented
    min_size : index
        the minimum number of vertices in a basin
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    depth_factor : float
        factor to determine whether to merge two neighboring watershed catchment
        basins -- they are merged if the Euclidean distance between their basin
        seeds is less than this fraction of the maximum Euclidean distance
        between points having minimum and maximum depths
    depth_ratio : float
        the minimum fraction of depth for a neighboring shallower
        watershed catchment basin (otherwise merged with the deeper basin)
    tolerance : float
        tolerance for detecting differences in depth between vertices
    regrow : bool
        regrow segments from watershed seeds?
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    segments : list of integers
        region numbers for all vertices
    seed_indices : list of integers
        list of indices to seed vertices

    Examples
    --------
    >>> # Perform watershed segmentation on the deeper portions of a surface:
    >>> import numpy as np
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.guts.segment import watershed
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> depth_file = fetch_data(urls['left_travel_depth'])
    >>> points, indices, f1, faces, depths, f2, npoints, f3 = read_vtk(depth_file,
    ...     return_first=True, return_array=True)
    >>> indices = np.where(depths > 0.01)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> min_size = 50
    >>> depth_factor = 0.25
    >>> depth_ratio = 0.1
    >>> tolerance = 0.01
    >>> regrow = True
    >>> background_value = -1
    >>> verbose = False
    >>> segments, seed_indices = watershed(depths, points,
    ...     indices, neighbor_lists, min_size, depth_factor, depth_ratio,
    ...     tolerance, regrow, background_value, verbose)
    >>> len(np.unique(segments))
    202
    >>> len_segments = []
    >>> for useg in np.unique(segments):
    ...     len_segments.append(len(np.where(segments == useg)[0]))
    >>> len_segments[0:10]
    [2976, 4092, 597, 1338, 1419, 1200, 1641, 220, 1423, 182]

    Write watershed segments and seeds to vtk file and view (skip test).
    Note: white spots indicate incomplete segmentation:

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> segments_array = np.array(segments)
    >>> segments_array[seed_indices] = 1.5 * np.max(segments_array)
    >>> rewrite_scalars(depth_file, 'watershed.vtk',
    ...                 segments_array, 'segments', [], -1) # doctest: +SKIP
    >>> plot_surfaces('watershed.vtk') # doctest: +SKIP

    """
    import numpy as np
    from time import time
    from mindboggle.guts.segment import extract_borders
    from mindboggle.guts.segment import segment_regions
    from mindboggle.guts.compute import point_distance

    # Make sure argument is a list
    if isinstance(indices, np.ndarray):
        indices = indices.tolist()

    if verbose:
        print('Segment {0} vertices by a surface watershed algorithm'.
              format(len(indices)))
        verbose2 = False
    else:
        verbose2 = False

    merge = True
    t0 = time()
    tiny = 0.000001

    use_depth_ratio = True

    #-------------------------------------------------------------------------
    # Find the borders of the given mesh vertices (indices):
    #-------------------------------------------------------------------------
    D = np.ones(len(depths))
    D[indices] = 2
    borders, foo1, foo2 = extract_borders(list(range(len(depths))), D,
        neighbor_lists, ignore_values=[], return_label_pairs=False)

    #-------------------------------------------------------------------------
    # Select deepest vertex as initial seed:
    #-------------------------------------------------------------------------
    index_deepest = indices[np.argmax(depths[indices])]
    seed_list = [index_deepest]
    basin_depths = []
    original_indices = indices[:]

    #-------------------------------------------------------------------------
    # Loop until all vertices have been segmented.
    # This limits the number of possible seeds:
    #-------------------------------------------------------------------------
    segments = background_value * np.ones(len(depths))
    seed_indices = []
    seed_points = []
    all_regions = []
    region = []
    counter = 0
    terminate = False
    while not terminate:

        # Add seeds to region:
        region.extend(seed_list)
        all_regions.extend(seed_list)

        # Remove seeds from vertices to segment:
        indices = list(frozenset(indices).difference(seed_list))
        if indices:

            # Identify neighbors of seeds:
            neighbors = []
            [neighbors.extend(neighbor_lists[x]) for x in seed_list]

            # Select neighbors that have not been previously selected
            # and are among the vertices to segment:
            old_seed_list = seed_list[:]
            seed_list = list(frozenset(neighbors).intersection(indices))
            seed_list = list(frozenset(seed_list).difference(all_regions))

            # For each vertex, select neighbors that are shallower:
            seed_neighbors = []
            for seed in old_seed_list:
                seed_neighbors.extend([x for x in neighbor_lists[seed]
                                       if depths[x] - tolerance <= depths[seed]])
            seed_list = list(frozenset(seed_list).intersection(seed_neighbors))

        else:
            seed_list = []

        # If there are no seeds remaining:
        if not len(seed_list):

            # If there is at least min_size points, assign counter to
            # segmented region, store index, and increment counter:
            if len(region) >= min_size:
                segments[region] = counter
                seed_indices.append(index_deepest)
                seed_points.append(points[index_deepest])
                counter += 1

                # Compute basin depth (max - min):
                Imax = region[np.argmax(depths[region])]
                Imin = region[np.argmin(depths[region])]
                max_depth = point_distance(points[Imax], [points[Imin]])[0]
                basin_depths.append(max_depth)

            # If vertices left to segment, re-initialize parameters:
            if indices:

                # Initialize new region/basin:
                region = []

                # Select deepest unsegmented vertex as new seed
                # if its rescaled depth is close to 1:
                index_deepest = indices[np.argmax(depths[indices])]
                seed_list = [index_deepest]

            # Termination criteria:
            if not len(indices):
                terminate = True

            # Display current number and size of region:
            if verbose2:
                print("    {0} vertices remain".format(len(indices)))

    if verbose:
        print('  ...Segmented {0} initial watershed regions ({1:.2f} seconds)'.
              format(counter, time() - t0))

    #-------------------------------------------------------------------------
    # Regrow from (deep) watershed seeds, stopping at borders:
    #-------------------------------------------------------------------------
    if regrow:

        if verbose:
            print('  Regrow segments from watershed seeds, '
                  'stopping at borders')
        indices = original_indices[:]
        segments = background_value * np.ones(len(depths))
        all_regions = []
        for iseed, seed_index in enumerate(seed_indices):
            seed_list = [seed_index]
            region = []
            terminate = False
            while not terminate:

                # Add seeds to region:
                region.extend(seed_list)
                all_regions.extend(seed_list)

                # Remove seeds from vertices to segment:
                indices = list(frozenset(indices).difference(seed_list))
                if indices:

                    # Identify neighbors of seeds:
                    neighbors = []
                    [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                    # Select neighbors that have not been previously selected
                    # and are among the vertices to segment:
                    old_seed_list = seed_list[:]
                    seed_list = list(frozenset(neighbors).intersection(indices))
                    seed_list = list(frozenset(seed_list).difference(all_regions))

                    # For each vertex, select neighbors that are shallower:
                    seed_neighbors = []
                    for seed in old_seed_list:
                        seed_neighbors.extend([x for x in neighbor_lists[seed]
                            if depths[x] - tolerance <= depths[seed]])
                    seed_list = list(frozenset(seed_list).intersection(seed_neighbors))

                    # Remove seed list if it contains a border vertex:
                    if seed_list:
                        if list(frozenset(seed_list).intersection(borders)):
                            seed_list = []
                else:
                    seed_list = []

                # Terminate growth for this seed if the seed_list is empty:
                if not len(seed_list):
                    terminate = True

                    # If there is at least min_size points, store index:
                    if len(region) >= min_size:
                        segments[region] = iseed

                    # Display current number and size of region:
                    if verbose2:
                        print("    {0} vertices remain".format(len(indices)))

        #---------------------------------------------------------------------
        # Continue growth until there are no more vertices to segment:
        #---------------------------------------------------------------------
        # Note: As long as keep_seeding=False, the segment values in `segments`
        # are equal to the order of the `basin_depths` and `seed_points` below.
        seed_lists = [[i for i,x in enumerate(segments) if x==s]
                      for s in np.unique(segments) if s != background_value]
        segments = segment_regions(indices, neighbor_lists, 1, seed_lists,
                                   False, False, [], [], [], '',
                                   background_value, False)

        if verbose:
            print('  ...Regrew {0} watershed regions from seeds '
                  '({1:.2f} seconds)'.format(iseed+1, time() - t0))

    #-------------------------------------------------------------------------
    # Merge watershed catchment basins:
    #-------------------------------------------------------------------------
    if merge:

        # Extract segments pairs at borders between watershed basins:
        if verbose:
            print('  Merge watershed catchment basins with deeper '
                  'neighboring basins')
            if verbose2:
                print('    Extract basin borders')
        foo1, foo2, pairs = extract_borders(original_indices, segments,
                                            neighbor_lists,
                                            ignore_values=[background_value],
                                            return_label_pairs=True)
        # Sort basin depths (descending order) -- return segment indices:
        Isort = np.argsort(basin_depths).tolist()
        Isort.reverse()

        # Find neighboring basins to each of the sorted basins:
        if verbose2:
            print("    Find neighboring basins")
        basin_pairs = []
        for index in Isort:
            index_neighbors = [int(list(frozenset(x).difference([index]))[0])
                               for x in pairs if index in x]
            if index_neighbors:

                # Store neighbors whose depth is less than a fraction of the
                # basin's depth and farther away than half the basin's depth:
                if use_depth_ratio:
                    index_neighbors = [[x, index] for x in index_neighbors
                        if basin_depths[x] / (basin_depths[index]+tiny) < depth_ratio
                        if point_distance(seed_points[x], [seed_points[index]])[0] >
                          depth_factor * max([basin_depths[x], basin_depths[index]])]
                # Store neighbors farther away than half the basin's depth:
                else:
                    index_neighbors = [[x, index] for x in index_neighbors
                        if point_distance(seed_points[x], [seed_points[index]])[0] >
                        depth_factor * max([basin_depths[x], basin_depths[index]])]
                if index_neighbors:
                    basin_pairs.extend(index_neighbors)

        # Merge shallow watershed catchment basins:
        if basin_pairs:
            if verbose2:
                print('    Merge basins with deeper neighboring basins')
            for basin_pair in basin_pairs:
                segments[np.where(segments == basin_pair[0])] = basin_pair[1]

        # Renumber segments so they are sequential:
        renumber_segments = segments.copy()
        segment_numbers = [int(x) for x in np.unique(segments)
                           if x != background_value]
        for i_segment, n_segment in enumerate(segment_numbers):
            segment = [i for i,x in enumerate(segments) if x == n_segment]
            renumber_segments[segment] = i_segment
        segments = renumber_segments

        # Print statement:
        if verbose:
            print('  ...Merged segments to form {0} watershed regions '
                  '({1:.2f} seconds)'.format(i_segment + 1, time() - t0))

    return segments.tolist(), seed_indices


def select_largest(points, faces, exclude_labels=[-1], areas=None,
                   reindex=True, background_value=-1, verbose=False):
    """
    Select the largest segment (connected set of indices) in a mesh.

    In case a surface patch is fragmented, we select the largest fragment,
    remove extraneous triangular faces, and reindex indices.

    Parameters
    ----------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh
    exclude_labels : list of integers
        background values to exclude
    areas : numpy array or list of floats (or None)
        surface area scalar values for all vertices
    reindex : bool
        reindex indices in faces?
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh

    Examples
    --------
    >>> # Two surface patches:
    >>> import numpy as np
    >>> from mindboggle.mio.vtks import read_scalars, read_vtk
    >>> from mindboggle.guts.mesh import keep_faces
    >>> from mindboggle.guts.segment import select_largest
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> exclude_labels = [-1]
    >>> points, indices, f1, faces, labels, f2,f3,f4 = read_vtk(label_file,
    ...                                                         True, True)
    >>> I28 = [i for i,x in enumerate(labels) if x==1028] # superior frontal
    >>> I20 = [i for i,x in enumerate(labels) if x==1020] # pars triangularis
    >>> I28.extend(I20)
    >>> faces = keep_faces(faces, I28)
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> reindex = True
    >>> background_value = -1
    >>> verbose = False
    >>> points2, faces2 = select_largest(points, faces, exclude_labels, areas,
    ...                                  reindex, background_value, verbose)
    >>> points2[0]
    [-1.4894376993179321, 53.260337829589844, 56.523414611816406]
    >>> points2[1]
    [-2.2537832260131836, 53.045711517333984, 56.23670959472656]
    >>> points2[2]
    [-13.091879844665527, 56.41604232788086, 59.330955505371094]
    >>> faces2[0]
    [7704, 7306, 7703]
    >>> faces2[1]
    [7694, 7704, 7705]
    >>> faces2[2]
    [7703, 8119, 7704]

    Write two surfaces to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import write_vtk # doctest: +SKIP
    >>> scalars = np.zeros(np.shape(labels)) # doctest: +SKIP
    >>> scalars[I28] = 1 # doctest: +SKIP
    >>> vtk_file = 'select_largest_two_labels.vtk' # doctest: +SKIP
    >>> write_vtk(vtk_file, points, indices, [], faces) # doctest: +SKIP
    >>> plot_surfaces(vtk_file) # doctest: +SKIP

    Write larger surface to vtk file and view (skip test):

    >>> vtk_file = 'select_largest.vtk' # doctest: +SKIP
    >>> write_vtk(vtk_file, points2, list(range(len(points2))), [], faces2) # doctest: +SKIP
    >>> plot_surfaces(vtk_file) # doctest: +SKIP

    """
    import numpy as np

    from mindboggle.guts.mesh import find_neighbors, keep_faces, \
        reindex_faces_points
    from mindboggle.guts.segment import segment_regions

    # Areas:
    use_area = False
    if isinstance(areas, np.ndarray) and np.shape(areas):
        use_area = True
    elif isinstance(areas, list) and len(areas):
        areas = np.array(areas)
        use_area = True

    # Check to see if there are enough points:
    min_npoints = 4
    npoints = len(points)
    if npoints < min_npoints or len(faces) < min_npoints:
        if verbose:
            print("The input size {0} ({1} faces) should be much larger "
                  "than {2}". format(npoints, len(faces), min_npoints))
        return None
    else:

        #---------------------------------------------------------------------
        # Segment the indices into connected sets of indices:
        #---------------------------------------------------------------------
        # Construct neighbor lists:
        neighbor_lists = find_neighbors(faces, npoints)

        # Determine the unique indices that make up the faces:
        indices = np.unique(np.ravel(faces))

        # Segment:
        segments = segment_regions(indices, neighbor_lists, 1, [], False,
                                   False, [], [], [], '', background_value,
                                   verbose)

        #---------------------------------------------------------------------
        # Select the largest segment (connected set of indices):
        #---------------------------------------------------------------------
        unique_segments = [x for x in np.unique(segments)
                           if x not in exclude_labels]
        if len(unique_segments) > 1:
            select_indices = []
            max_segment_area = 0
            for segment_number in unique_segments:
                segment_indices = [i for i,x in enumerate(segments)
                                   if x == segment_number]
                if use_area:
                    segment_area = np.sum(areas[segment_indices])
                else:
                    segment_area = len(segment_indices)
                if segment_area > max_segment_area:
                    select_indices = segment_indices
                    max_segment_area = segment_area

                # Print message:
                if verbose:
                    if use_area:
                        print('Segment {0}: {1} vertices ({2:.2f} area)'.
                              format(int(segment_number),
                                     len(segment_indices),
                                     segment_area))
                    else:
                        print('Segment {0}: {1} vertices'.
                              format(int(segment_number),
                                     len(segment_indices)))
            if verbose:
                print('Largest of {0} segments: {1:.2f}'.
                      format(len(unique_segments), max_segment_area))

            #-----------------------------------------------------------------
            # Renumber faces for the selected indices:
            #-----------------------------------------------------------------
            faces = keep_faces(faces, select_indices)
            if faces:
                #-------------------------------------------------------------
                # Reindex indices in faces:
                #-------------------------------------------------------------
                if reindex:
                    faces, points, o1 = reindex_faces_points(faces, points)
                else:
                    points = np.array(points)
                    points = points[select_indices].tolist()
                return points, faces
            else:
                return None
        else:
            return points, faces


def extract_borders(indices, labels, neighbor_lists,
                    ignore_values=[], return_label_pairs=False):
    """
    Detect the label borders in a collection of vertices such as a region.

    Label borders are the set of all vertices
    whose neighbors do not share the same label.

    Parameters
    ----------
    indices : list of integers
        indices to (a subset of) vertices
    labels : numpy array of integers
        label numbers for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    ignore_values : list of integers
        integers to ignore (e.g., background)
    return_label_pairs : bool
        return label pairs?

    Returns
    -------
    border_indices : list of integers
        indices to label boundary vertices
    border_label_tuples : list of lists of sorted pairs of integers
        sorted label pairs
    unique_border_label_tuples : list of sorted pairs of integers
        unique, sorted label pairs

    Examples
    --------
    >>> # Small example:
    >>> from mindboggle.guts.segment import extract_borders
    >>> indices = [0,1,2,4,5,8,9]
    >>> labels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, -1, -1]
    >>> neighbor_lists = [[1,2,3], [1,2], [2,3], [2], [4,7], [3,2,3],
    ...                   [2,3,7,8], [2,6,7], [3,4,8], [7], [7,8], [3,2,3]]
    >>> border_indices, border_label_tuples, ubLT = extract_borders(indices,
    ...     labels, neighbor_lists, [], True)
    >>> border_indices
    [0, 1, 2, 4, 5, 8]
    >>> border_label_tuples
    [[20, 30, 40], [20, 30], [30, 40], [50, 80], [30, 40], [40, 50, 90]]
    >>> ubLT
    [[20, 30, 40], [20, 30], [30, 40], [50, 80], [40, 50, 90]]

    Real example -- extract sulcus label boundaries:

    >>> import numpy as np
    >>> from mindboggle.guts.mesh import find_neighbors
    >>> from mindboggle.guts.segment import extract_borders
    >>> from mindboggle.mio.vtks import read_vtk
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> f1,f2,f3, faces, labels, f4, npoints, f5 = read_vtk(label_file,
    ...                                                     True, True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> indices_borders, label_pairs, f1 = extract_borders(list(range(npoints)),
    ...     labels, neighbor_lists)
    >>> indices_borders[0:10]
    [115, 116, 120, 121, 125, 126, 130, 131, 281, 286]

    Write borders on surfaces to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars
    >>> IDs = -1 * np.ones(npoints) # doctest: +SKIP
    >>> IDs[indices_borders] = 1 # doctest: +SKIP
    >>> rewrite_scalars(label_file, 'extract_borders.vtk',
    ...                 IDs, 'borders') # doctest: +SKIP
    >>> plot_surfaces('extract_borders.vtk') # doctest: +SKIP

    Write just the borders to vtk file and view (skip test):

    >>> rewrite_scalars(label_file, 'extract_borders_no_background.vtk',
    ...                 IDs, 'borders', IDs) # doctest: +SKIP
    >>> plot_surfaces('extract_borders_no_background.vtk') # doctest: +SKIP

    """
    import numpy as np

    # Make sure arguments are numpy arrays:
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)

    # Construct an array of labels corresponding to the neighbor lists:
    L = np.array([list(set(labels[lst])) for lst in neighbor_lists])

    # Find indices to sets of two labels:
    border_indices = [indices[y] for y in
                      [i for i,x in enumerate(L[indices])
                       if len(set(x)) >= 2]]

    if return_label_pairs:
        border_label_tuples = [np.sort(L[indices[j]]).tolist() for j in
                                [i for i,x in enumerate(L[indices])
                                 if len(set(x)) >= 2]]
    else:
        border_label_tuples = []

    if ignore_values:
        Ikeep = [i for i,x in enumerate(border_label_tuples)
                 if not len(set(x).intersection(ignore_values))]
        border_indices = [x for i,x in enumerate(border_indices)
                            if i in Ikeep]
        if return_label_pairs:
            border_label_tuples = [x for i,x in enumerate(border_label_tuples)
                                    if i in Ikeep]

    if return_label_pairs:
        unique_border_label_tuples = []
        for pair in border_label_tuples:
            if pair not in unique_border_label_tuples:
                unique_border_label_tuples.append(pair)
    else:
        unique_border_label_tuples = []

    return border_indices, border_label_tuples, unique_border_label_tuples


def extract_borders_2nd_surface(labels_file, values_file='',
                                output_file='', background_value=-1):
    """
    Extract borders (between labels) on a surface.
    Option: extract border values on a second surface.

    Parameters
    ----------
    labels_file : string
        file name for surface mesh with labels
    values_file : string
        file name for surface mesh with values to extract along borders
    output_file : string
        full path to output file
    background_value : integer or float
        background value

    Returns
    -------
    border_file : string
        full path to output file
    border_values : numpy array
        values for all vertices (background for vertices off label borders)
    indices_borders : list of integers
        indices to borders vertices

    Examples
    --------
    >>> # Extract depth values along label borders in sulci (mask):
    >>> import numpy as np
    >>> from mindboggle.guts.segment import extract_borders_2nd_surface
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> values_file = fetch_data(urls['left_travel_depth'])
    >>> background_value = -1
    >>> output_file = 'extract_borders_2nd_surface.vtk'
    >>> border_file, values, I = extract_borders_2nd_surface(label_file,
    ...     values_file, output_file, background_value)
    >>> print(np.array_str(np.unique(values)[0:8],
    ...       precision=5, suppress_small=True))
    [-1.       0.       0.00012  0.00023  0.00032  0.00044  0.00047  0.00051]
    >>> I[0:10]
    [115, 116, 120, 121, 125, 126, 130, 131, 281, 286]

    Write depth values on label borders to vtk file and view (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> plot_surfaces(border_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    from mindboggle.mio.vtks import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.guts.mesh import find_neighbors
    from mindboggle.guts.segment import extract_borders

    # Load labeled surface file
    points, indices, lines, faces, labels, scalar_names, npoints, \
        input_vtk = read_vtk(labels_file, return_first=True,
                             return_array=True)

    # Detect borders
    neighbor_lists = find_neighbors(faces, npoints)
    indices_borders, foo1, foo2 = extract_borders(list(range(npoints)),
                                                  labels, neighbor_lists)

    # Filter values with label borders
    border_values = background_value * np.ones(npoints)
    if values_file:
        values, name = read_scalars(values_file, return_first=True, return_array=True)
        border_values[indices_borders] = values[indices_borders]
    else:
        border_values[indices_borders] = 1

    # Write out label boundary vtk file
    if output_file:
        border_file = output_file
    else:
        border_file = os.path.join(os.getcwd(),
                                   'borders_' + os.path.basename(labels_file))
    rewrite_scalars(labels_file, border_file, border_values, \
                    'label_borders_in_mask', [], background_value)

    if not os.path.exists(border_file):
        raise IOError(border_file + " not found")

    return border_file, border_values, indices_borders


def combine_2labels_in_2volumes(file1, file2, label1=3, label2=2,
                                output_file=''):
    """
    Combine the voxels of one label from two files and overlay them
    on the second label combined from the two files.

    An example application is to combine two segmentation volumes,
    such as from FreeSurfer and ANTs, to obtain a single cortex=2
    and noncortex=3 segmentation file, by taking the union of cortex voxels
    from the segmentations, the union of the noncortex voxels from the
    segmentations, and overwriting intersecting cortex and noncortex voxels
    with noncortex (3) labels.

    ANTs tends to include more cortical gray matter at the periphery of
    the brain than Freesurfer, and FreeSurfer tends to include more white
    matter that extends deep into gyral folds than ANTs, so this function
    attempts to remedy their differences by overlaying ANTs cortical gray
    with FreeSurfer white matter.

    Example preprocessing steps ::

      1. Run Freesurfer and antsCorticalThickness.sh on T1-weighted image.
      2. Convert FreeSurfer volume labels (e.g., wmparc.mgz or aparc+aseg.mgz)
         to cortex (2) and noncortex (3) segments using relabel_volume()
         function [refer to LABELS.rst or FreeSurferColorLUT labels file].
      3. Convert ANTs Atropos-segmented volume (tmpBrainSegmentation.nii.gz)
         to cortex and noncortex segments, by converting 1-labels to 0 and
         4-labels to 3 with the relabel_volume() function
         (the latter is to include deep-gray matter with noncortical tissues).

    Parameters
    ----------
    file1 : string
        name of nibabel-readable file with labels in {label1, label2}
    file2 : string
        name of nibabel-readable file with labels in {label1, label2}
    label1 : integer
        source label number
    label2 : integer
        target label number (to be overwritten by label1 where they intersect)
    output_file : string
        output file name

    Returns
    -------
    output_file : string
        name of combined segmented nibabel-readable file
        with labels in {label1, label2}

    Examples
    --------
    >>> # Example following documentation above:
    >>> import os
    >>> from mindboggle.guts.segment import combine_2labels_in_2volumes
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> file1 = fetch_data(urls['freesurfer_segmentation'])
    >>> file2 = fetch_data(urls['ants_segmentation'])
    >>> label1 = 3
    >>> label2 = 2
    >>> output_file = 'combine_2labels_in_2volumes.nii.gz'
    >>> output_file = combine_2labels_in_2volumes(file1, file2, label1,
    ...                                           label2, output_file)

    View nifti file (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes(output_file) # doctest: +SKIP

    """
    import os
    import numpy as np
    import nibabel as nb

    #-------------------------------------------------------------------------
    # Load labeled image volume and extract data as 1-D array:
    #-------------------------------------------------------------------------
    vol1 = nb.load(file1)
    vol2 = nb.load(file2)
    data1 = vol1.get_data().ravel()
    data2 = vol2.get_data().ravel()
    xfm = vol1.get_affine()
    #-------------------------------------------------------------------------
    # Indices to voxels with label1 or label2 in two files:
    #-------------------------------------------------------------------------
    label1 = int(label1)
    label2 = int(label2)
    Ilabel1_vol1 = np.where(data1 == label1)[0]
    Ilabel1_vol2 = np.where(data2 == label1)[0]
    Ilabel2_vol1 = np.where(data1 == label2)[0]
    Ilabel2_vol2 = np.where(data2 == label2)[0]
    Ilabel1 = Ilabel1_vol1.tolist() + Ilabel1_vol2.tolist()
    Ilabel2 = Ilabel2_vol1.tolist() + Ilabel2_vol2.tolist()
    #-------------------------------------------------------------------------
    # Assign new labels and reshape to original dimensions:
    #-------------------------------------------------------------------------
    new_data = data1.copy()
    new_data[Ilabel2] = label2
    new_data[Ilabel1] = label1
    new_data = np.reshape(new_data, vol1.shape)
    #-------------------------------------------------------------------------
    # Save relabeled file:
    #-------------------------------------------------------------------------
    if not output_file:
        output_file = os.path.join(os.getcwd(),
                                   'combined_segmentations.nii.gz')
    img = nb.Nifti1Image(new_data, xfm)
    img.to_filename(output_file)

    if not os.path.exists(output_file):
        raise IOError(output_file + " not found")

    return output_file


def split_brain(image_file, label_file, left_labels, right_labels):
    """
    Split a brain using left/right labels.

    Parameters
    ----------
    image_file : string
        nibabel-readable image volume
    label_file : string
        nibabel-readable labeled brain image volume
    left_labels : list of integers
        left label numbers
    right_labels : list of integers
        right label numbers

    Returns
    -------
    left_brain : string
        name of nibabel-readable image volume with left half of brain image
    right_brain : string
        name of nibabel-readable image volume with right half of brain image

    Examples
    --------
    >>> from mindboggle.guts.segment import split_brain
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> image_file = fetch_data(urls['freesurfer_segmentation'])
    >>> label_file = fetch_data(urls['freesurfer_labels'])
    >>> dkt = DKTprotocol()
    >>> left_labels = dkt.left_cerebrum_numbers
    >>> right_labels = dkt.right_cerebrum_numbers
    >>> left_brain, right_brain = split_brain(image_file, label_file,
    ...                                       left_labels, right_labels)

    Write results to nifti file and view (skip test):

    >>> from mindboggle.mio.plots import plot_volumes # doctest: +SKIP
    >>> plot_volumes([left_brain, right_brain]) # doctest: +SKIP

    """
    import os
    import numpy as np
    import nibabel as nb

    from mindboggle.guts.relabel import keep_volume_labels

    left_brain = os.path.join(os.getcwd(),
                              'left_' + os.path.basename(image_file))
    right_brain = os.path.join(os.getcwd(),
                               'right_' + os.path.basename(image_file))
    #-------------------------------------------------------------------------
    # Split brain labels by masking with left or right labels:
    #-------------------------------------------------------------------------
    left_brain = keep_volume_labels(label_file, left_labels, left_brain)
    right_brain = keep_volume_labels(label_file, right_labels, right_brain)
    #-------------------------------------------------------------------------
    # Load labeled image volumes and extract data as 1-D array:
    #-------------------------------------------------------------------------
    vol = nb.load(image_file)
    volL = nb.load(left_brain)
    volR = nb.load(right_brain)
    data = vol.get_data().ravel()
    dataL = volL.get_data().ravel()
    dataR = volR.get_data().ravel()
    dataL[np.where(dataL != 0)[0]] = 1
    dataR[np.where(dataR != 0)[0]] = 1
    xfm = vol.get_affine()
    #-------------------------------------------------------------------------
    # Split brain image by masking with left or right labels:
    #-------------------------------------------------------------------------
    left_data = data * dataL
    right_data = data * dataR
    #-------------------------------------------------------------------------
    # Save relabeled file:
    #-------------------------------------------------------------------------
    left_data = np.reshape(left_data, volL.shape)
    right_data = np.reshape(right_data, volR.shape)
    img1 = nb.Nifti1Image(left_data, xfm)
    img1.to_filename(left_brain)
    img2 = nb.Nifti1Image(right_data, xfm)
    img2.to_filename(right_brain)

    if not os.path.exists(right_brain) or not os.path.exists(left_brain):
        raise IOError(right_brain + " or " + left_brain + "not found")

    return left_brain, right_brain


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules