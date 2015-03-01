#!/usr/bin/env python
"""
Functions for segmenting a surface mesh.

Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Propagate labels to segment surface into contiguous regions
#------------------------------------------------------------------------------
def propagate(points, faces, region, seeds, labels,
              max_iters=500, tol=0.001, sigma=10):
    """
    Propagate labels to segment surface into contiguous regions,
    starting from seed vertices.

    Imports : mindboggle.labels.rebound

    Parameters
    ----------
    points : array (or list) of lists of three integers
        coordinates for all vertices
    faces : list of lists of three integers
        indices to three vertices per face (indices start from zero)
    region : numpy array of integers
        values > -1 indicate inclusion in a region for all vertices
    seeds : numpy array of integers
        seed numbers for all vertices (default -1 for not a seed)
    labels : numpy array of integers
        label numbers for all vertices, with -1s for unlabeled vertices
    max_iters : integer
        maximum number of iterations to run graph-based learning algorithm
    tol: float
        threshold to assess convergence of the algorithm
    sigma: float
        gaussian kernel parameter

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> # Propagate labels between label boundary segments in a single fold:
    >>> import os
    >>> import numpy as np
    >>> import mindboggle.labels.rebound as rb
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.label import extract_borders, propagate
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> indices_boundaries, label_pairs, foo = extract_borders(range(npoints),
    >>>     labels, neighbor_lists)
    >>> label_pair_lists = sulcus_boundaries()
    >>> # Select a single fold
    >>> fold_ID = 2
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>> fold_array = -1 * np.ones(npoints)
    >>> fold_array[indices_fold] = 1
    >>> # Extract the boundary for this fold
    >>> indices_boundaries, label_pairs, foo = extract_borders(indices_fold,
    >>>     labels, neighbor_lists)
    >>> # Select boundary segments in the sulcus labeling protocol
    >>> seeds = -1 * np.ones(npoints)
    >>> for ilist,label_pair_list in enumerate(label_pair_lists):
    >>>     I = [x for i,x in enumerate(indices_boundaries) if np.sort(label_pairs[i]).tolist() in label_pair_list]
    >>>     seeds[I] = ilist
    >>> #
    >>> segments = propagate(points, faces, fold_array, seeds, labels)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(labels_file, 'test_propagate.vtk',
    >>>                 segments, 'segments', segments)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_propagate.vtk')

    """
    import numpy as np
    from mindboggle.utils.mesh import remove_faces
    import mindboggle.utils.kernels as kernels
    import mindboggle.labels.rebound as rb

    indices_region = [i for i,x in enumerate(region) if x > -1]
    if indices_region and faces:
        local_indices_region = -1 * np.ones(region.size)
        local_indices_region[indices_region] = range(len(indices_region))

        # Make sure arguments are numpy arrays
        if not isinstance(seeds, np.ndarray):
            seeds = np.array(seeds)
        if not isinstance(labels, np.ndarray):
            labels = np.array(labels)
        if not isinstance(points, np.ndarray):
            points = np.array(points)

        if points.size:

            n_sets = len(np.unique([x for x in seeds if x > -1]))
            if n_sets == 1:
                print('Segment {0} vertices from 1 set of seed vertices'.
                      format(len(indices_region)))
            else:
                print('Segment {0} vertices from {1} sets of seed vertices'.
                      format(len(indices_region), n_sets))

            # Set up rebound Bounds class instance
            B = rb.Bounds()
            B.Faces = remove_faces(faces, indices_region)
            if B.Faces:
                B.Indices = local_indices_region
                B.Points = points[indices_region]
                B.Labels = labels[indices_region]
                B.seed_labels = seeds[indices_region]
                B.num_points = len(B.Points)

                # Propagate seed IDs from seeds
                B.graph_based_learning(method='propagate_labels', realign=False,
                                       kernel=kernels.rbf_kernel, sigma=sigma,
                                       max_iters=max_iters, tol=tol, vis=False)
            else:
                print("  No faces")

            # Assign maximum probability seed IDs to each point of region
            max_prob_labels = B.assign_max_prob_label()

            # Return segment IDs in original vertex array
            segments = -1 * np.ones(len(points))
            segments[indices_region] = max_prob_labels
    else:
        segments = -1 * np.ones(len(points))

    return segments

#------------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions by seed growing
#------------------------------------------------------------------------------
def segment(vertices_to_segment, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False,
            spread_within_labels=False, labels=[], label_lists=[], values=[]):
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
    keep_seeding : Boolean
        grow from new seeds even after all seed lists have fully grown
    spread_within_labels : Boolean
        grow seeds only by vertices with labels in the seed labels?
    labels : list of integers (required only if spread_within_labels)
        label numbers for all vertices, with -1s for unlabeled vertices
    label_lists : list of lists of integers (required only if spread_within_labels)
        List of unique labels for each seed list to grow into
        (If empty, set to unique labels for each seed list)
    values : list of floats (default empty)
        values for all vertices for use in preferentially directed segmentation
        (segment in direction of lower values)

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> # Segment deep regions with or without seeds:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.segment import segment
    >>> from mindboggle.labels.label import extract_borders
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> vertices_to_segment = np.where(depths > 0.50)[0].tolist()  # higher to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> # Example 1: without seed lists
    >>> folds = segment(vertices_to_segment, neighbor_lists)
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_segment.vtk', folds, 'folds', folds)
    >>> plot_vtk('test_segment.vtk')
    >>> #
    >>> # Example 2: with seed lists
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> label_pair_lists = sulcus_boundaries()
    >>> label_lists = [np.unique(np.ravel(x)) for x in label_pair_lists]
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file)
    >>> indices_boundaries, label_pairs, foo = extract_borders(vertices_to_segment,
    >>>     labels, neighbor_lists)
    >>> seed_lists = []
    >>> for label_pair_list in label_pair_lists:
    >>>     seed_lists.append([x for i,x in enumerate(indices_boundaries) if np.sort(label_pairs[i]).tolist() in label_pair_list])
    >>> #
    >>> sulci = segment(vertices_to_segment, neighbor_lists, 1,
    >>>                 seed_lists, True, True, labels, label_lists, values=[])
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_segment_seeds.vtk', sulci, 'sulci', sulci)
    >>> plot_vtk('test_segment_seeds.vtk')

    """
    import numpy as np

    verbose = False

    # If seed_lists is empty, select first vertex from vertices_to_segment
    # (single vertex selection does not affect result -- see below*)
    if seed_lists:
        select_single_seed = False
        if len(seed_lists) == 1:
            print('    Segment {0} vertices using seed vertices'.
                  format(len(vertices_to_segment)))
        else:
            print('    Segment {0} vertices from {1} sets of seed vertices'.
                  format(len(vertices_to_segment), len(seed_lists)))
    else:
        select_single_seed = True
        seed_lists = [[vertices_to_segment[0]]]
        print('    Segment {0} vertices (first vertex as initial seed)'.
              format(len(vertices_to_segment)))

    # Initialize variables, including the list of vertex indices for each region,
    # vertex indices for all regions, and Boolean list indicating which regions
    # are not fully grown, max. region ID, number of segments, etc.
    segments = -1 * np.ones(len(neighbor_lists))
    region_lists = [[] for x in seed_lists]
    all_regions = []
    fully_grown = [False for x in seed_lists]
    new_segment_index = 0
    counter = 0

    # If label_lists empty, set to unique labels for each seed list
    if spread_within_labels:
        if not len(label_lists):
            label_lists = []
            for seed_list in seed_lists:
                seed_labels = np.unique([labels[x] for x in seed_list])
                label_lists.append(seed_labels)

    # Loop until all of the seed lists have grown to their full extent
    while not all(fully_grown):
        # Loop through seed lists over and over again
        for ilist, seed_list in enumerate(seed_lists):
            # If seed list empty
            if not len(seed_list):
                fully_grown[ilist] = True
                # If seed list not fully grown
            if not fully_grown[ilist]:

                # Add seeds to region
                region_lists[ilist].extend(seed_list)
                all_regions.extend(seed_list)

                # Remove seeds from vertices to segment
                vertices_to_segment = list(frozenset(vertices_to_segment).
                difference(seed_list))

                if vertices_to_segment:

                    # Find neighbors of each seed with lower values than the seed
                    if values:
                        neighbors = []
                        for seed in seed_list:
                            seed_neighbors = neighbor_lists[seed]
                            seed_neighbors = [x for x in seed_neighbors
                                              if values[x] <= values[seed]]
                            if seed_neighbors:
                                neighbors.extend(seed_neighbors)
                    else:
                        # Find neighbors of seeds
                        neighbors = []
                        [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                    # Select neighbors that have not been previously selected
                    # and are among the vertices to segment
                    seed_list = list(frozenset(neighbors).intersection(vertices_to_segment))
                    seed_list = list(frozenset(seed_list).difference(all_regions))

                else:
                    seed_list = []

                # If there are seeds remaining
                if seed_list:

                    # Select neighbors with the same labels
                    # as the initial seed labels
                    if spread_within_labels:
                        seed_list = [x for x in seed_list
                                     if labels[x] in label_lists[ilist]]

                    # Continue growing seed list
                    seed_lists[ilist] = seed_list

                # If there are no seeds remaining
                else:
                    # Stop growing seed list (see exception below)
                    fully_grown[ilist] = True

                    # If the region size is large enough
                    size_region = len(region_lists[ilist])
                    if size_region >= min_region_size:

                        # Assign ID to segmented region and increment ID
                        if select_single_seed:
                            new_segment_index = counter
                            counter += 1
                        else:
                            new_segment_index = ilist
                        segments[region_lists[ilist]] = new_segment_index

                        # Display current number and size of region
                        if verbose and size_region > 1:
                            if len(seed_lists) == 1:
                                if vertices_to_segment:
                                    print("      {0} vertices remain".
                                          format(len(vertices_to_segment)))
                            else:
                                print("      Region {0}: {1} vertices ({2} remain)".
                                      format(int(new_segment_index), size_region,
                                             len(vertices_to_segment)))

                    # If selecting a single seed, continue growing
                    # if there are more vertices to segment
                    if select_single_seed:
                        if len(vertices_to_segment) >= min_region_size:
                            fully_grown[0] = False
                            seed_lists[0] = [vertices_to_segment[0]]
                            region_lists[0] = []

    # Keep growing from new seeds even after all seed lists have fully grown
    if keep_seeding and len(vertices_to_segment) >= min_region_size:
        print('    Keep seeding to segment {0} remaining vertices'.
              format(len(vertices_to_segment)))

        # Select first unsegmented vertex as new seed
        seed_list = [vertices_to_segment[0]]

        # Loop until the seed list has grown to its full extent
        new_segment_index = ilist + 1
        region = []
        while len(vertices_to_segment) >= min_region_size:

            # Add seeds to region
            region.extend(seed_list)
            all_regions.extend(seed_list)

            # Remove seeds from vertices to segment
            vertices_to_segment = list(frozenset(vertices_to_segment).
            difference(seed_list))
            if vertices_to_segment:

                # Identify neighbors of seeds
                neighbors = []
                [neighbors.extend(neighbor_lists[x]) for x in seed_list]

                # Select neighbors that have not been previously selected
                # and are among the vertices to segment
                seed_list = list(frozenset(vertices_to_segment).intersection(neighbors))
                seed_list = list(frozenset(seed_list).difference(all_regions))
            else:
                seed_list = []

            # If there are no seeds remaining
            if not len(seed_list):

                # If the region size is large enough
                size_region = len(region)
                if size_region >= min_region_size:

                    # Assign ID to segmented region and increment ID
                    segments[region] = new_segment_index
                    new_segment_index += 1

                    # Display current number and size of region
                    if verbose and size_region > 1:
                        print("      {0} vertices remain".
                              format(len(vertices_to_segment)))

                # Select first unsegmented vertex as new seed
                if len(vertices_to_segment) >= min_region_size:
                    seed_list = [vertices_to_segment[0]]
                    region = []

    return segments

#------------------------------------------------------------------------------
# Fill boundaries on a surface mesh to segment vertices into contiguous regions
#------------------------------------------------------------------------------
def segment_by_filling_boundaries(regions, neighbor_lists):
    """
    Fill boundaries (contours) on a surface mesh
    to segment vertices into contiguous regions.

    Steps ::
        1. Extract region boundaries (assumed to be closed contours)
        2. Segment boundaries into separate, contiguous boundaries
        3. For each boundary
            4. Find the neighbors to either side of the boundary
            5. Segment the neighbors into exterior and interior sets of neighbors
            6. Find the interior (smaller) sets of neighbors
            7. Fill the contours formed by the interior neighbors

    Parameters
    ----------
    regions : numpy array of integers
        region numbers for all vertices (default -1)
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> # Segment folds by extracting their boundaries and filling them in separately:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.segment import segment_by_filling_boundaries
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> regions = -1 * np.ones(npoints)
    >>> regions[depths > 0.50] = 1
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> folds = segment_by_filling_boundaries(regions, neighbor_lists)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_segment_by_filling_boundaries.vtk', folds, 'folds', folds)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_segment_by_filling_boundaries.vtk')

    """
    import numpy as np
    from mindboggle.labels.label import extract_borders, segment

    include_boundary = False

    # Make sure arguments are numpy arrays
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    print('Segment vertices using region boundaries')

    # Extract region boundaries (assumed to be closed contours)
    print('  Extract region boundaries (assumed to be closed contours)')
    indices_boundaries, label_pairs, foo = extract_borders(range(len(regions)),
                                                           regions, neighbor_lists)
    # Extract background
    indices_background = list(frozenset(range(len(regions))).
    difference(indices_boundaries))

    # Segment boundaries into separate, contiguous boundaries
    print('  Segment boundaries into separate, contiguous boundaries')
    boundaries = segment(indices_boundaries, neighbor_lists, 1)

    # For each boundary
    unique_boundaries = [x for x in np.unique(boundaries) if x > -1]
    segments = -1 * np.ones(len(regions))
    for boundary_number in unique_boundaries:

        print('  Boundary {0} of {1}:'.format(int(boundary_number),
                                              len(unique_boundaries)))
        boundary_indices = [i for i,x in enumerate(boundaries)
                            if x == boundary_number]
        # Find the neighbors to either side of the boundary
        indices_neighbors = []
        [indices_neighbors.extend(neighbor_lists[i]) for i in boundary_indices]
        #indices_neighbors2 = indices_neighbors[:]
        #[indices_neighbors2.extend(neighbor_lists[i]) for i in indices_neighbors]
        indices_neighbors = list(frozenset(indices_neighbors).
        difference(indices_boundaries))

        # Segment the neighbors into exterior and interior sets of neighbors
        print('    Segment the neighbors into exterior and interior sets of neighbors')
        neighbors = segment(indices_neighbors, neighbor_lists, 1)

        # Find the interior (smaller) sets of neighbors
        print('    Find the interior (smaller) sets of neighbors')
        seed_lists = []
        unique_neighbors = [x for x in np.unique(neighbors) if x > -1]
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
        print('    Fill the contour formed by the interior neighbors')
        vertices_to_segment = list(frozenset(indices_background).
        difference(indices_boundaries))
        segment_region = segment(vertices_to_segment, neighbor_lists, 1, [seed_list])

        if include_boundary:
            segment_region[boundary_indices] = 1

        segments[segment_region > -1] = boundary_number

    return segments

#------------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions by a watershed algorithm
#------------------------------------------------------------------------------
def watershed(depths, indices, neighbor_lists, depth_ratio=0.1, tolerance=0.01):
    """
    Segment vertices of surface into contiguous "watershed basin" regions
    by seed growing from an iterative selection of the deepest vertices.

    Note: This function, when used alone, has the same drawback as segment() --
    the order of seed selection influences the result for multiple
    seeds within a connected set of vertices (region).
    To ameliorate this bias, we run the shrink_segments() function
    on the segments returned by watershed(), which shrinks segments in
    regions with multiple segments, and use these fractional segments
    as seeds for the propagate() function, which is slower and insensitive
    to depth, but is not biased by seed order.

    Parameters
    ----------
    depths : numpy array of floats
        depth values for all vertices (default -1)
    indices : list of integers
        indices to mesh vertices to be segmented
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    depth_ratio : float
        minimum fraction of depth for a neighboring shallower watershed
        catchment basin (otherwise merged with the deeper basin)
    tolerance : float
        tolerance for detecting differences in depth between vertices

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices (default -1)

    Examples
    --------
    >>> # Perform watershed segmentation on the deeper portions of a surface:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.segment import watershed
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> indices = np.where(depths > 0.11)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> depth_ratio = 0.1
    >>> tolerance = 0.01
    >>> #
    >>> segments = watershed(depths, indices, neighbor_lists,
    >>>                      depth_ratio, tolerance)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_watershed.vtk',
    >>>                 segments, 'segments', segments)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_watershed.vtk')

    """
    import numpy as np
    from time import time
    from mindboggle.labels.label import extract_borders

    print('Segment {0} vertices by a surface watershed algorithm'.
          format(len(indices)))
    verbose = False
    merge = True
    t0 = time()

    # Select deepest vertex as initial seed
    seed_list = [indices[np.argmax(depths[indices])]]
    #minima = [seed_list[0]]
    max_depths = [depths[seed_list[0]]]
    basin_depths = []
    original_indices = indices[:]

    # Loop until all vertices segmented
    segments = -1 * np.ones(len(depths))
    all_regions = []
    region = []
    counter = 0
    while len(indices):

        if verbose:
            print('  Segment {0} vertices (deepest vertex as initial seed)'.
                  format(len(indices)))

        # Add seeds to region
        region.extend(seed_list)
        all_regions.extend(seed_list)

        # Remove seeds from vertices to segment
        indices = list(frozenset(indices).difference(seed_list))
        if indices:

            # Identify neighbors of seeds
            neighbors = []
            [neighbors.extend(neighbor_lists[x]) for x in seed_list]

            # Select neighbors that have not been previously selected
            # and are among the vertices to segment
            old_seed_list = seed_list[:]
            seed_list = list(frozenset(neighbors).intersection(indices))
            seed_list = list(frozenset(seed_list).difference(all_regions))

            # For each vertex, select neighbors that are shallower
            seed_neighbors = []
            for seed in old_seed_list:
                seed_neighbors.extend([x for x in neighbor_lists[seed]
                                       if depths[x] - tolerance <= depths[seed]])
            seed_list = list(frozenset(seed_list).intersection(seed_neighbors))

        else:
            seed_list = []

        # If there are no seeds remaining
        if not len(seed_list):

            # Assign ID to segmented region and increment ID
            segments[region] = counter

            # Compute basin depth (max - min)
            basin_depths.append(max_depths[-1] - np.min(depths[region]))

            # If vertices left to segment, re-initialize parameters
            if indices:
                # Select deepest unsegmented vertex as new seed
                seed_list = [indices[np.argmax(depths[indices])]]
                # Initialize new region/basin, maximum depth, and counter
                region = []
                max_depths.append(depths[seed_list[0]])
                #minima.append(seed_list[0])
                counter += 1

            # Display current number and size of region
            if verbose:
                print("    {0} vertices remain".format(len(indices)))
    print('  ...Segmented {0} regions ({1:.2f} seconds)'.
          format(counter, time() - t0))

    # Merge watershed catchment basins
    if merge:

        # Extract boundaries between watershed catchment basins
        print('  Merge watershed catchment basins with deeper neighboring basins')
        print('    Extract region boundaries')
        foo1, foo2, pairs = extract_borders(original_indices, segments,
                                              neighbor_lists, ignore_indices=[-1])

        # Sort basin depths (descending order) -- return segment indices
        Isort = np.argsort(basin_depths).tolist()
        Isort.reverse()

        # Find neighboring basins to each of the sorted basins
        print("    Find neighbors whose depth is less than a fraction of the basin's depth")
        basin_pairs = []
        for index in Isort:
            index_neighbors = [int(list(frozenset(x).difference([index]))[0])
                               for x in pairs if index in x]
            # Store neighbors whose depth is less than a fraction of the basin's depth
            tiny = 0.000001
            index_neighbors = [[x, index] for x in index_neighbors
                               if basin_depths[x] / (basin_depths[index] + tiny) < depth_ratio]
            if index_neighbors:
                basin_pairs.extend(index_neighbors)

        # Merge shallow watershed catchment basins
        if basin_pairs:
            print('    Merge basins with deeper neighboring basins')
            for basin_pair in basin_pairs:
                segments[np.where(segments == basin_pair[0])] = basin_pair[1]

        # Print statement
        n_segments = len([x for x in np.unique(segments) if x > -1])
        print('  ...Segmented and merged {0} watershed regions ({1:.2f} seconds)'.
              format(n_segments, time() - t0))

    return segments

#------------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions by seed growing
#------------------------------------------------------------------------------
def shrink_segments(regions, segments, depths, remove_fraction=0.75,
                    only_multiple_segments=False):
    """
    Shrink segments in a segmented surface mesh by a fraction of its maximum
    depth, for all segments or for segments in regions with multiple segments.

    The segment() and watershed() functions, when used alone, are influenced
    by the order of seed selection for multiple seeds within a connected
    set of vertices (region, such as a fold).
    To ameliorate this bias, we run this function on the segments returned by
    segment() or watershed() to shrink segments in regions with multiple
    segments (only_multiple_segments=True).
    The output can be used as seeds for the propagate() function,
    which is not biased by seed order.

    Parameters
    ----------
    regions : list or array of integers
        region IDs for all vertices, indicating inclusion in a region (default -1)
    segments : numpy array of integers
        segment IDs for all vertices, indicating inclusion in a segment (default -1)
    depths : numpy array of floats
        depth values for all vertices (default -1)
    remove_fraction : float
        remove fraction of the previously segmented depth
        for regions of connected vertices with multiple segments
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
    >>> from mindboggle.labels.segment import watershed, shrink_segments
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno',
    >>>                                      'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> indices = np.where(depths > 0.11)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> segments = watershed(depths, indices, neighbor_lists,
    >>>                      depth_ratio=0.1, tolerance=0.01)
    >>> #
    >>> shrunken_segments = shrink_segments(folds, segments, depths,
    >>>     remove_fraction=0.25, only_multiple_segments=True)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'test_shrink_segments.vtk',
    >>>     shrunken_segments, 'shrunken_segments', shrunken_segments)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_shrink_segments.vtk')

    """
    import numpy as np

    print('Shrink segments')

    shrunken_segments = -1 * np.ones(len(depths))

    # Make sure arguments are numpy arrays
    if not isinstance(segments, np.ndarray):
        segments = np.array(segments)
    if not isinstance(depths, np.ndarray):
        depths = np.array(depths)

    # Shrink only segments in regions with multiple segments
    if only_multiple_segments:
        print('  Remove {0:.2f} depth from each segment for regions with '
              'multiple segments'.format(remove_fraction))

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
        print('  Remove {0:.2f} of each segment by depth'.format(remove_fraction))
        unique_segments = [x for x in np.unique(segments) if x > -1]
        for n_segment in unique_segments:
            indices_segment = [i for i,x in enumerate(segments) if x == n_segment]
            depth_threshold = remove_fraction * np.max(depths[indices_segment])
            indices_segment = [x for x in indices_segment
                               if depths[x] > depth_threshold]
            shrunken_segments[indices_segment] = n_segment

    return shrunken_segments

#------------------------------------------------------------------------------
# Find vertices with highest values within a fraction of the surface
#------------------------------------------------------------------------------
def extract_high_values(values, areas, fraction):
    """
    Find the highest-valued vertices in a surface whose collective area
    is a given fraction of the total surface area of the mesh.

    Example: extract half of a surface with the highest depth values

    Parameters
    ----------
    values : numpy array of floats
        values for all vertices
    areas : numpy array of floats
        surface area values for all vertices
    fraction : float
        fraction of surface mesh to extract

    Returns
    -------
    area_values : array of integers
        an integer for every mesh vertex: 1 if extracted, -1 if not

    Examples
    --------
    >>> # Extract the highest depth values constituting half of the surface area:
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> from mindboggle.labels.segment import extract_high_values
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(path, 'arno', 'measures', 'lh.pial.area.vtk')
    >>> values, name = read_scalars(vtk_file, return_first=True, return_array=True)
    >>> areas, name = read_scalars(area_file, return_first=True, return_array=True)
    >>> fraction = 0.50
    >>> #
    >>> area_values = extract_high_values(values, areas, fraction)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(vtk_file, 'test_extract_high_values.vtk',
    >>>                 area_values, 'area values', area_values)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_high_values.vtk')

    """
    import numpy as np

    print("  Extract the highest-valued surface vertices ({0} of surface area)".
          format(fraction))

    # Make sure arguments are numpy arrays
    if not isinstance(values, np.ndarray):
        values = np.array(values)
    if not isinstance(areas, np.ndarray):
        areas = np.array(areas)

    # Load depth and surface area values from VTK files
    indices_asc = np.argsort(values)
    indices_des = indices_asc[::-1]

    total_area = np.sum(areas)
    fraction_area = fraction * total_area
    area_values = -1 * np.ones(len(areas))

    # Start with fraction_area of the vertices
    start = np.round(fraction * len(values))
    area_values[indices_des[0:start]] = 1
    sum_area = np.sum(areas[indices_des[0:start]])

    # If these initial vertices cover less than fraction_area,
    # add vertices until the remaining vertices' area exceeds fraction_area
    if sum_area <= fraction_area:
        for index in indices_des[start::]:
            area_values[index] = 1
            sum_area += areas[index]
            if sum_area >= fraction_area:
                break
    # Otherwise, if these initial vertices cover more than fraction_area,
    # remove vertices until the remaining vertices' area is less than fraction_area
    else:
        start = np.round((1-fraction) * len(values))
        for index in indices_asc[start::]:
            area_values[index] = -1
            sum_area += areas[index]
            if sum_area <= fraction_area:
                break

    return area_values
