#!/usr/bin/env python
"""
Functions for segmenting a surface mesh.

Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-----------------------------------------------------------------------------
# Propagate ring_labels to segment surface into contiguous regions
#-----------------------------------------------------------------------------
def propagate(points, faces, region, seeds, labels,
              max_iters=500, tol=0.001, sigma=10, background_value=-1):
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
    region : list (or array) of integers
        values > background_value indicate inclusion in a region for all vertices
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

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Propagate labels between label boundary segments in a single fold:
    >>> import os
    >>> import numpy as np
    >>> import mindboggle.labels.rebound as rb
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.labels.labels import extract_borders
    >>> from mindboggle.utils.segment import propagate
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    >>> from mindboggle.labels.protocol import dkt_protocol
    >>> protocol = 'DKT25'
    >>> sulcus_names, sulcus_label_pair_lists, unique_sulcus_label_pairs,
    ...     label_names, label_numbers, cortex_names, cortex_numbers,
    ...     noncortex_names, noncortex_numbers = dkt_protocol(protocol)
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> folds, name = read_scalars(folds_file, return_first=True, return_array=True)
    >>> faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file,
    >>>     return_first=True, return_array=True)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> indices_borders, label_pairs, foo = extract_borders(range(npoints),
    >>>     labels, neighbor_lists)
    >>> background_value = -1
    >>> # Select a single fold
    >>> fold_ID = 2
    >>> indices_fold = [i for i,x in enumerate(folds) if x == fold_ID]
    >>> fold_array = background_value * np.ones(npoints)
    >>> fold_array[indices_fold] = 1
    >>> # Extract the boundary for this fold
    >>> indices_borders, label_pairs, foo = extract_borders(indices_fold,
    >>>     labels, neighbor_lists)
    >>> # Select boundary segments in the sulcus labeling protocol
    >>> seeds = background_value * np.ones(npoints)
    >>> for ilist,label_pair_list in enumerate(sulcus_label_pair_lists):
    >>>     I = [x for i,x in enumerate(indices_borders)
    >>>          if np.sort(label_pairs[i]).tolist() in label_pair_list]
    >>>     seeds[I] = ilist
    >>> #
    >>> segments = propagate(points, faces, fold_array, seeds, labels)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(labels_file, 'propagate.vtk',
    >>>                 segments, 'segments', segments)
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> plot_surfaces('propagate.vtk')

    """
    import numpy as np
    from mindboggle.utils.mesh import remove_faces
    import mindboggle.utils.kernels as kernels
    import mindboggle.labels.rebound as rb

    # Make sure arguments are numpy arrays
    if not isinstance(seeds, np.ndarray):
        seeds = np.array(seeds)
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)
    if not isinstance(points, np.ndarray):
        points = np.array(points)

    indices_region = [i for i,x in enumerate(region) if x != background_value]
    if indices_region and faces:
        local_indices_region = background_value * np.ones(labels.shape)
        local_indices_region[indices_region] = range(len(indices_region))

        if points.size:

            n_sets = len(np.unique([x for x in seeds if x != background_value]))
            if n_sets == 1:
                print('Segment {0} vertices from 1 set of seed vertices'.
                      format(len(indices_region)))
            else:
                print('Segment {0} vertices from {1} sets of seed vertices'.
                      format(len(indices_region), n_sets))

            # Remove faces whose three vertices are not among specified indices:
            refaces = remove_faces(faces, indices_region)

            # Set up rebound Bounds class instance
            B = rb.Bounds()
            if refaces:
                B.Faces = np.array(refaces)
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
            segments = background_value * np.ones(len(points))
            segments[indices_region] = max_prob_labels
    else:
        segments = background_value * np.ones(len(points))

    return segments

#-----------------------------------------------------------------------------
# Segment vertices of surface into contiguous regions by seed growing
#-----------------------------------------------------------------------------
def segment(vertices_to_segment, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False, spread_within_labels=False,
            labels=[], label_lists=[], values=[], max_steps='', verbose=False,
            background_value=-1):
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
        label numbers for all vertices
    label_lists : list of lists of integers (required only if spread_within_labels)
        List of unique labels for each seed list to grow into
        (If empty, set to unique labels for each seed list)
    values : list of floats (default empty)
        values for all vertices for use in preferentially directed segmentation
        (segment in direction of lower values)
    max_steps : integer (or empty string for infinity)
        maximum number of segmentation steps to take for each seed list
    verbose : Boolean
        print out text?
    background_value : integer
        background value

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Segment deep regions with or without seeds:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.utils.segment import segment
    >>> from mindboggle.labels.labels import extract_borders
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> vertices_to_segment = np.where(depths > 0.50)[0].tolist()  # higher to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> # Example 1: without seed lists
    >>> folds = segment(vertices_to_segment, neighbor_lists)
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'segment.vtk', folds, 'folds', folds)
    >>> plot_surfaces('segment.vtk')
    >>> #
    >>> # Example 2: with seed lists
    >>> from mindboggle.labels.protocol import dkt_protocol
    >>> protocol = 'DKT25'
    >>> sulcus_names, sulcus_label_pair_lists, unique_sulcus_label_pairs,
    ...     label_names, label_numbers, cortex_names, cortex_numbers,
    ...     noncortex_names, noncortex_numbers = dkt_protocol(protocol)
    >>> label_lists = [np.unique(np.ravel(x)) for x in sulcus_label_pair_lists]
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, labels, name, input_vtk = read_vtk(labels_file)
    >>> indices_borders, label_pairs, foo = extract_borders(vertices_to_segment,
    >>>     labels, neighbor_lists, ignore_values=[], return_label_pairs=True)
    >>> seed_lists = []
    >>> for label_pair_list in sulcus_label_pair_lists:
    >>>     seed_lists.append([x for i,x in enumerate(indices_borders) if np.sort(label_pairs[i]).tolist() in label_pair_list])
    >>> #
    >>> sulci = segment(vertices_to_segment, neighbor_lists, 1,
    >>>                 seed_lists, True, True, labels, label_lists, values=[])
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'segment_seeds.vtk', sulci, 'sulci', sulci)
    >>> plot_surfaces('segment_seeds.vtk')

    """
    import numpy as np

    verbose = False

    # Make sure arguments are lists:
    if isinstance(vertices_to_segment, np.ndarray):
        vertices_to_segment = vertices_to_segment.tolist()
    if isinstance(labels, np.ndarray):
        labels = labels.tolist()
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

#-----------------------------------------------------------------------------
# Fill borders on a surface mesh to segment vertices into contiguous regions
#-----------------------------------------------------------------------------
def segment_by_filling_borders(regions, neighbor_lists, background_value=-1):
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
    background_value : integer
        background value

    Returns
    -------
    segments : numpy array of integers
        region numbers for all vertices

    Examples
    --------
    >>> # Segment folds by extracting their borders and filling them in separately:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.utils.segment import segment_by_filling_borders
    >>> from mindboggle.utils.io_vtk import read_vtk, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> background_value = -1
    >>> regions = background_value * np.ones(npoints)
    >>> regions[depths > 0.50] = 1
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> #
    >>> folds = segment_by_filling_borders(regions, neighbor_lists, background_value)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'segment_by_filling_borders.vtk', folds, 'folds', folds)
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> plot_surfaces('segment_by_filling_borders.vtk')

    """
    import numpy as np
    from mindboggle.labels.labels import extract_borders
    from mindboggle.utils.segment import segment

    include_boundary = False

    # Make sure arguments are numpy arrays
    if not isinstance(regions, np.ndarray):
        regions = np.array(regions)

    print('Segment vertices using region borders')

    # Extract region borders (assumed to be closed contours)
    print('  Extract region borders (assumed to be closed contours)')
    indices_borders, foo1, foo2 = extract_borders(range(len(regions)),
                                        regions, neighbor_lists)
    # Extract background
    indices_background = list(frozenset(range(len(regions))).
    difference(indices_borders))

    # Segment borders into separate, contiguous borders
    print('  Segment borders into separate, contiguous borders')
    borders = segment(indices_borders, neighbor_lists, 1)

    # For each boundary
    unique_borders = [x for x in np.unique(borders) if x != background_value]
    segments = background_value * np.ones(len(regions))
    for boundary_number in unique_borders:

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
        print('    Segment the neighbors into exterior and interior sets of neighbors')
        neighbors = segment(indices_neighbors, neighbor_lists, 1)

        # Find the interior (smaller) sets of neighbors
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
        print('    Fill the contour formed by the interior neighbors')
        vertices_to_segment = list(frozenset(indices_background).
        difference(indices_borders))
        segment_region = segment(vertices_to_segment, neighbor_lists, 1, [seed_list])

        if include_boundary:
            segment_region[border_indices] = 1

        segments[segment_region != background_value] = boundary_number

    return segments

def segment_rings(region, seeds, neighbor_lists, step=1, background_value=-1):
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

    Returns
    -------
    segments : list of lists of integers
        indices to vertices for each concentric segment

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.labels.labels import extract_borders
    >>> from mindboggle.utils.segment import segment_rings
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> values_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> values, name = read_scalars(values_file, True, True)
    >>> vtk_file = os.path.join(path, 'arno', 'freesurfer', 'lh.pial.vtk')
    >>> neighbor_lists = find_neighbors_from_file(vtk_file)
    >>> background_value = -1
    >>> fold_file = os.path.join(path, 'arno', 'features', 'fold11.vtk')
    >>> fold, name = read_scalars(fold_file)
    >>> indices = [i for i,x in enumerate(fold) if x != background_value]
    >>> # Initialize seeds with the boundary of thresholded indices:
    >>> use_threshold = True
    >>> if use_threshold:
    >>>     # Threshold at the median depth or within maximum values in boundary:
    >>>     threshold = np.median(values[indices]) #+ np.std(values[indices])
    >>>     indices_high = [x for x in indices if values[x] >= threshold]
    >>>     # Make sure threshold is within the maximum values of the boundary:
    >>>     B = np.ones(len(values))
    >>>     B[indices] = 2
    >>>     borders, foo1, foo2 = extract_borders(range(len(B)), B, neighbor_lists)
    >>>     borders = [x for x in borders if values[x] != background_value]
    >>>     if list(frozenset(indices_high).intersection(borders)):
    >>>         threshold = np.max(values[borders]) + np.std(values[borders])
    >>>         indices_high = [x for x in indices if values[x] >= threshold]
    >>>     # Extract threshold boundary vertices as seeds:
    >>>     B = background_value * np.ones(len(values))
    >>>     B[indices_high] = 2
    >>>     seeds, foo1, foo2 = extract_borders(range(len(values)), B, neighbor_lists)
    >>> # Or initialize P with the maximum value point:
    >>> else:
    >>>     seeds = [indices[np.argmax(values[indices])]]
    >>>     indices_high = []
    >>> #
    >>> indices = list(frozenset(indices).difference(indices_high))
    >>> indices = list(frozenset(indices).difference(seeds))
    >>> segments = segment_rings(indices, seeds, neighbor_lists, step=1)
    >>> #
    >>> # View:
    >>> S = background_value * np.ones(len(values))
    >>> for i, segment in enumerate(segments):
    >>>     S[segment] = i
    >>> rewrite_scalars(vtk_file, 'segment_rings.vtk', S, 'segment_rings', fold)
    >>> plot_surfaces('segment_rings.vtk')
    >>> # Store:
    >>> #import pickle
    >>> #output_file = os.path.join(path, 'tests', 'segments_fold11.pkl')
    >>> #pickle.dump(segments, open(output_file, "wb" ))

    """
    from mindboggle.utils.segment import segment

    segments = []
    while seeds:

        # Segment step-wise starting from seeds and through the region:
        seeds_plus_new = segment(region, neighbor_lists, min_region_size=1,
                                 seed_lists=[seeds], keep_seeding=False,
                                 spread_within_labels=False, labels=[],
                                 label_lists=[], values=[], max_steps=step)
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
# Segment vertices of surface into contiguous regions by a watershed algorithm
#-----------------------------------------------------------------------------
def watershed(depths, points, indices, neighbor_lists, min_size=1,
              depth_factor=0.25, depth_ratio=0.1, tolerance=0.01, regrow=True,
              background_value=-1):
    """
    Segment vertices of a surface mesh into contiguous "watershed basins"
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
    regrow : Boolean
        regrow segments from watershed seeds?
    background_value : integer
        background value

    Returns
    -------
    segments : list of integers
        region numbers for all vertices
    seed_indices : list of integers
        list of indices to seed vertices

    Examples
    --------
    >>> # Perform watershed segmentation on the deeper portions of a surface:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.mesh import find_neighbors
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> from mindboggle.utils.segment import watershed, segment
    >>> from mindboggle.utils.io_vtk import read_vtk, read_scalars, rewrite_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> faces, lines, indices, points, npoints, depths, name, input_vtk = read_vtk(depth_file,
    >>>     return_first=True, return_array=True)
    >>> indices = np.where(depths > 0.01)[0]  # high to speed up
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> min_size = 50
    >>> depth_factor = 0.25
    >>> depth_ratio = 0.1
    >>> tolerance = 0.01
    >>> regrow = True
    >>> background_value = -1
    >>> #
    >>> segments, seed_indices = watershed(depths, points,
    >>>     indices, neighbor_lists, min_size, depth_factor, depth_ratio,
    >>>     tolerance, regrow, background_value)
    >>> #
    >>> # Write results to vtk file and view:
    >>> rewrite_scalars(depth_file, 'watershed.vtk',
    >>>                 segments, 'segments', segments)
    >>> plot_surfaces('watershed.vtk')
    >>> # View watershed seeds:
    >>> seeds = background_value * np.ones(len(depths))
    >>> for i, s in enumerate(seed_indices):
    >>>     seeds[s] = i
    >>> rewrite_scalars(depth_file, 'watershed_seeds.vtk',
    >>>                 seeds, 'seeds', seeds)
    >>> plot_surfaces('watershed_seeds.vtk')

    """
    import numpy as np
    from time import time
    from mindboggle.labels.labels import extract_borders
    from mindboggle.utils.segment import segment
    from mindboggle.utils.compute import point_distance

    # Make sure argument is a list
    if isinstance(indices, np.ndarray):
        indices.tolist()

    print('Segment {0} vertices by a surface watershed algorithm'.
          format(len(indices)))
    verbose = False
    merge = True
    t0 = time()
    tiny = 0.000001

    use_depth_ratio = True

    #-------------------------------------------------------------------------
    # Find the borders of the given mesh vertices (indices):
    #-------------------------------------------------------------------------
    D = np.ones(len(depths))
    D[indices] = 2
    borders, foo1, foo2 = extract_borders(range(len(depths)), D,
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
            if verbose:
                print("    {0} vertices remain".format(len(indices)))

    print('  ...Segmented {0} initial watershed regions ({1:.2f} seconds)'.
          format(counter, time() - t0))

    #-------------------------------------------------------------------------
    # Regrow from (deep) watershed seeds, stopping at borders:
    #-------------------------------------------------------------------------
    if regrow:

        print('  Regrow segments from watershed seeds, stopping at borders')
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
                    if verbose:
                        print("    {0} vertices remain".format(len(indices)))

        #---------------------------------------------------------------------
        # Continue growth until there are no more vertices to segment:
        #---------------------------------------------------------------------
        # Note: As long as keep_seeding=False, the segment values in `segments`
        # are equal to the order of the `basin_depths` and `seed_points` below.
        seed_lists = [[i for i,x in enumerate(segments) if x==s]
                      for s in np.unique(segments) if s != background_value]
        segments = segment(indices, neighbor_lists, min_region_size=1,
            seed_lists=seed_lists, keep_seeding=False, spread_within_labels=False,
            labels=[], label_lists=[], values=[], max_steps='', verbose=False)

        print('  ...Regrew {0} watershed regions from seeds ({1:.2f} seconds)'.
              format(iseed+1, time() - t0))

    #-------------------------------------------------------------------------
    # Merge watershed catchment basins:
    #-------------------------------------------------------------------------
    if merge:

        # Extract segments pairs at borders between watershed basins:
        print('  Merge watershed catchment basins with deeper neighboring basins')
        if verbose:
            print('    Extract basin borders')
        foo1, foo2, pairs = extract_borders(original_indices, segments,
                                            neighbor_lists,
                                            ignore_values=[background_value],
                                            return_label_pairs=True)
        # Sort basin depths (descending order) -- return segment indices:
        Isort = np.argsort(basin_depths).tolist()
        Isort.reverse()

        # Find neighboring basins to each of the sorted basins:
        if verbose:
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
            if verbose:
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
        print('  ...Merged segments to form {0} watershed regions ({1:.2f} seconds)'.
              format(i_segment + 1, time() - t0))

    return segments.tolist(), seed_indices


def select_largest(points, faces, exclude_labels=[-1], areas=None,
                   reindex=True):
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
    reindex : Boolean
        reindex indices in faces?

    Returns
    -------
    points : list of lists of 3 floats
        x,y,z coordinates for each vertex of the structure
    faces : list of lists of 3 integers
        3 indices to vertices that form a triangle on the mesh

    Examples
    --------
    >>> # Two surface patches:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    >>> from mindboggle.utils.mesh import remove_faces
    >>> from mindboggle.utils.segment import select_largest
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> exclude_labels = [-1]
    >>> faces, lines, indices, points, u1, labels, u2,u3 = read_vtk(label_file,
    >>>      return_first=True, return_array=True)
    >>> I28 = [i for i,x in enumerate(labels) if x==28] # superior frontal
    >>> I20 = [i for i,x in enumerate(labels) if x==20] # pars triangularis
    >>> I28.extend(I20)
    >>> faces = remove_faces(faces, I28)
    >>> areas, u1 = read_scalars(area_file, True, True)
    >>> reindex = True
    >>> #
    >>> points2, faces2 = select_largest(points, faces, exclude_labels, areas,
    >>>                                  reindex)
    >>> # View:
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> scalars = np.zeros(np.shape(labels))
    >>> scalars[I28] = 1
    >>> vtk_file = 'test_two_labels.vtk'
    >>> write_vtk(vtk_file, points, indices, lines, faces,
    >>>           scalars, scalar_names='scalars')
    >>> plot_surfaces(vtk_file)

    """
    import numpy as np

    from mindboggle.utils.mesh import find_neighbors, remove_faces, \
        reindex_faces_points
    from mindboggle.utils.segment import segment

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
        segments = segment(indices, neighbor_lists, min_region_size=1,
            seed_lists=[], keep_seeding=False, spread_within_labels=False,
            labels=[], label_lists=[], values=[], max_steps='', verbose=False)

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
                if use_area:
                    print('Segment {0}: {1} vertices ({2:.2f} area)'.
                          format(int(segment_number), len(segment_indices),
                                 segment_area))
                else:
                    print('Segment {0}: {1} vertices'.
                          format(int(segment_number), len(segment_indices)))
            print('Largest of {0} segments: {1:.2f}'.
                  format(len(unique_segments), max_segment_area))

            #-----------------------------------------------------------------
            # Renumber faces for the selected indices:
            #-----------------------------------------------------------------
            faces = remove_faces(faces, select_indices)
            if faces:
                #-------------------------------------------------------------
                # Reindex indices in faces:
                #-------------------------------------------------------------
                if reindex:
                    faces, points = reindex_faces_points(faces, points)
                else:
                    points = np.array(points)
                    points = points[select_indices].tolist()
                return points, faces
            else:
                return None
        else:
            return points, faces
