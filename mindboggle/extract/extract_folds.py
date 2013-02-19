#!/usr/bin/env python
"""
Functions to extract folds, sulci, or identify sulci from folds.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

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
    depth_file : str
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
    >>> from mindboggle.utils.io_vtk import read_faces_points, rewrite_vtk
    >>> from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    >>> from mindboggle.extract.extract_folds import extract_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'measures', 'lh.pial.depth.vtk')
    >>> faces, points, npoints = read_faces_points(depth_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>>
    >>> folds, n_folds = extract_folds(depth_file, neighbor_lists, 50, True)
    >>>
    >>> # Write results to vtk file and view with mayavi2:
    >>> folds = folds.tolist()
    >>> rewrite_vtk(depth_file, 'test_extract_folds.vtk', folds, 'folds', folds)
    >>> os.system('mayavi2 -m Surface -d test_extract_folds.vtk &')

    """
    import numpy as np
    from time import time
    from scipy.ndimage.filters import gaussian_filter1d
    from mindboggle.utils.io_vtk import read_vtk
    from mindboggle.utils.mesh_operations import segment, fill_holes, \
        watershed, shrink_segments
    if not len(neighbor_lists):
        from mindboggle.utils.mesh_operations import find_neighbors

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
    if not len(neighbor_lists):
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
    """
    >>> # Plot histogram:
    >>> a,b,c = hist(depths, bins=nbins)
    """

    #---------------------------------------------------------------------------
    # Anticipating that there will be a rapidly decreasing distribution
    # of low depth values (on the outer surface) with a long tail of higher
    # depth values (in the folds), smooth the bin values (Gaussian), convolve
    # to compute slopes, and find the depth for the first bin with slope = 0.
    #---------------------------------------------------------------------------
    bins_smooth = gaussian_filter1d(bins.tolist(), 5)
    """
    >>> # Plot smoothed histogram:
    >>> plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')
    """
    window = [-1, 0, 1]
    bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
    ibin = np.where(bin_slopes == 0)[0]
    if len(ibin):
        depth_threshold = bin_edges[ibin[0]]
    else:
        depth_threshold = np.median(depths)

    #---------------------------------------------------------------------------
    # Find the deepest vertices
    #---------------------------------------------------------------------------
    indices_deep = [i for i,x in enumerate(depths) if x >= depth_threshold]
    if len(indices_deep):

        #-----------------------------------------------------------------------
        # Segment deep vertices as an initial set of folds
        #-----------------------------------------------------------------------
        print("  Segment vertices deeper than {0:.2f} as folds".format(depth_threshold))
        t1 = time()
        folds = segment(indices_deep, neighbor_lists)
        # Slightly slower alternative -- fill boundaries:
        #regions = -1 * np.ones(len(points))
        #regions[indices_deep] = 1
        #folds = fill_boundaries(regions, neighbor_lists)
        print('    ...Segmented folds ({0:.2f} seconds)'.format(time() - t1))
        """
        >>> # Display resulting initial folds:
        >>> rewrite_vtk(depth_file, 'test_folds1.vtk',
        >>>                      [folds.tolist()], 'folds', folds)
        >>> os.system('mayavi2 -m Surface -d test_folds1.vtk &')
        """

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
        do_fill_holes = True
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
# Extract sulci
#===============================================================================
def extract_sulci(surface_vtk, folds, labels, neighbor_lists, label_pair_lists,
                  min_boundary=1, sulcus_names=[]):
    """
    Identify sulci from folds in a brain surface according to a labeling
    protocol that includes a list of label pairs defining each sulcus.

    Definitions ::
        ``fold``: a group of connected, deep vertices
        ``label list``: list of labels used to define a single sulcus
        ``label pair list``: list of pairs of labels, where each pair
                             defines a boundary between two labeled regions;
                             no two label pair lists share a label pair
        ``sulcus ID``: ID to uniquely identify a sulcus,
                       as index to a sulcus label pair list

    Steps for each fold::
        1. Remove fold if it has fewer than two labels.
        2. Remove fold if its labels do not contain a sulcus label pair.
        3. Find vertices with labels that are in only one of the fold's
           label boundary pairs. Assign the vertices the sulcus with the
           label pair if they are connected to the label boundary for that pair.
        4. If there are remaining vertices, segment into sets of vertices
           connected to label boundaries, and assign a sulcus ID to each segment.

    Parameters
    ----------
    surface_vtk : string
        file name for surface mesh vtk (from which to extract points and faces)
    labels : list of integers
        labels for all vertices
    folds : list or array of integers
        fold IDs for all vertices
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    min_boundary : integer
        minimum number of vertices for a sulcus label boundary segment
    sulcus_names : list of strings (optional)
        names of sulci

    Returns
    -------
    sulci : array of integers
        sulcus numbers for all vertices, with -1s for non-sulcus vertices
    n_sulci : integers
        number of sulci

    Examples
    --------
    >>> import os
    >>> from time import time
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_vtk, write_vtk
    >>> from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    >>> from mindboggle.extract.extract_folds import extract_sulci
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>>
    >>> # Load labels, folds, neighbor lists, and sulcus names and label pairs
    >>> folds_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                      'features', 'lh.folds.vtk')
    >>> labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                            'labels', 'lh.labels.DKT25.manual.vtk')
    >>> faces, lines, indices, points, npoints, folds, name = read_vtk(folds_file)
    >>> faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file)
    >>> neighbor_lists = find_neighbors(faces, npoints)
    >>> sulcus_names_file = os.path.join(data_path, 'info', 'sulcus_names.txt')
    >>> fid = open(sulcus_names_file, 'r')
    >>> sulcus_names = fid.readlines()
    >>> sulcus_names = [x.strip('\n') for x in sulcus_names]
    >>> label_pair_lists = sulcus_boundaries()
    >>> min_boundary = 10
    >>>
    >>> # Extract sulci
    >>> sulci, n_sulci = extract_sulci(labels_file, folds, labels,
    >>>                                neighbor_lists, label_pair_lists,
    >>>                                min_boundary, sulcus_names)
    >>>
    >>> # Finally, write points, faces and sulci to a new vtk file
    >>> #rewrite_vtk(labels_file, 'test_extract_sulci.vtk',
    >>> #                     [sulci.tolist()], ['sulci'], sulci.tolist())
    >>> indices = [i for i,x in enumerate(sulci) if x > -1]
    >>> write_vtk('test_extract_sulci.vtk', points, indices, lines,
    >>>    inside_faces(faces, indices), [sulci.tolist()], ['sulci'])
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')

    """
    from time import time
    import numpy as np
    from mindboggle.utils.io_vtk import read_faces_points
    from mindboggle.utils.mesh_operations import detect_boundaries, propagate, segment

    #---------------------------------------------------------------------------
    # Prepare data
    #---------------------------------------------------------------------------
    # Array of sulcus IDs for fold vertices, initialized as -1.
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in the label list, or vertices having only one label,
    # their sulcus IDs will remain -1.
    sulci = -1 * np.ones(len(neighbor_lists))

    # Prepare list of sulcus label lists (labels within a sulcus)
    label_lists = []
    for row in label_pair_lists:
        label_lists.append(list(np.unique(np.asarray(row))))

    # Prepare list of all unique sorted label pairs in the labeling protocol
    protocol_pairs = []
    [protocol_pairs.append(np.unique(x).tolist()) for lst in label_pair_lists
     for x in lst if np.unique(x).tolist() not in protocol_pairs]

    # Load points, faces
    faces, points, npoints = read_faces_points(surface_vtk)

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    fold_numbers = [int(x) for x in np.unique(folds) if x > -1]
    n_folds = len(fold_numbers)
    print("Extract sulci from {0} folds...".format(n_folds))
    t0 = time()
    for n_fold in fold_numbers:
        fold = [i for i,x in enumerate(folds) if x == n_fold]
        len_fold = len(fold)
        # List the labels in this fold (greater than zero)
        fold_labels = [labels[x] for x in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels) if x > 0]

        #-----------------------------------------------------------------------
        # NO MATCH -- fold has fewer than two labels
        #-----------------------------------------------------------------------
        if len(unique_fold_labels) < 2:
            # Ignore: sulci already initialized with -1 values
            if len(unique_fold_labels) == 0:
                print("  Fold {0} ({1} vertices): NO MATCH -- fold has no labels".
                      format(n_fold, len_fold))
            else:
                print("  Fold {0} ({1} vertices): "
                  "NO MATCH -- fold has only one label ({2})".
                  format(n_fold, len_fold, unique_fold_labels[0]))
            # Ignore: sulci already initialized with -1 values

        else:
            # Find all label boundary pairs within the fold
            indices_fold_pairs, fold_pairs, unique_fold_pairs = detect_boundaries(
                fold, labels, neighbor_lists)

            # Find fold label pairs in the protocol (pairs are already sorted)
            fold_pairs_in_protocol = [x for x in unique_fold_pairs
                                      if x in protocol_pairs]

            print("  Fold {0} labels: {1} ({2} vertices)".format(n_fold,
                  ', '.join([str(x) for x in unique_fold_labels]), len_fold))
            print("  Fold {0} label pairs in protocol: {1}".format(n_fold,
                  ', '.join([str(x) for x in fold_pairs_in_protocol])))

            #-------------------------------------------------------------------
            # NO MATCH -- fold has no sulcus label pair
            #-------------------------------------------------------------------
            if not len(fold_pairs_in_protocol):
                print("  Fold {0}: NO MATCH -- fold has no sulcus label pair".
                      format(n_fold, len_fold))

            #-------------------------------------------------------------------
            # Possible matches
            #-------------------------------------------------------------------
            else:
                # Labels in the protocol (includes repeats across label pairs)
                labels_in_pairs = [x for lst in fold_pairs_in_protocol for x in lst]

                # Labels that appear in one or more than one sulcus label boundary
                unique_labels = []
                nonunique_labels = []
                for label in np.unique(labels_in_pairs):
                    if len([x for x in labels_in_pairs if x == label]) == 1:
                        unique_labels.append(label)
                    else:
                        nonunique_labels.append(label)

                #---------------------------------------------------------------
                # Vertices whose labels are in only one sulcus label pair
                #---------------------------------------------------------------
                # Find vertices with a label that is in only one of the fold's
                # label pairs (the other label in the pair can exist
                # in other pairs). Assign the vertices the sulcus with the label
                # pair if they are connected to the label boundary for that pair.
                #---------------------------------------------------------------
                if len(unique_labels):

                    for pair in fold_pairs_in_protocol:
                        # If one or both labels in label pair is/are unique
                        unique_labels_in_pair = [x for x in pair if x in unique_labels]
                        n_unique = len(unique_labels_in_pair)
                        if n_unique:

                            ID = [i for i,x in enumerate(label_pair_lists)
                                  if pair in x][0]

                            # Construct seeds from label boundary vertices
                            # (fold_pairs and pair already sorted)
                            indices_pair = [x for i,x in enumerate(indices_fold_pairs)
                                            if fold_pairs[i] == pair]

                            # Identify vertices with unique label(s) in pair
                            indices_unique_labels = [fold[i]
                                                     for i,x in enumerate(fold_labels)
                                                     if x in unique_labels_in_pair]

                            # Propagate from seeds to labels in label pair
                            sulci2 = segment(indices_unique_labels, neighbor_lists,
                                             min_region_size=1,
                                             seed_lists=[indices_pair],
                                             keep_seeding=False,
                                             spread_within_labels=True,
                                             labels=labels)
                            sulci[sulci2 > -1] = ID

                            # Print statement
                            if n_unique == 1:
                                ps1 = '1 label'
                            else:
                                ps1 = 'Both labels'
                            if len(sulcus_names):
                                ps2 = sulcus_names[ID]
                            else:
                                ps2 = ''
                            print("    {0} unique to one fold pair: {1} {2}".
                                  format(ps1, ps2, unique_labels_in_pair))

                #---------------------------------------------------------------
                # Vertex labels shared by multiple label pairs
                #---------------------------------------------------------------
                # Propagate labels from label boundaries to vertices with labels
                # that are shared by multiple label pairs in the fold.
                #---------------------------------------------------------------
                if len(nonunique_labels):
                    # For each label shared by different label pairs
                    for label in nonunique_labels:
                        # Print statement
                        print("    Propagate sulcus label boundaries with label {0}".
                              format(int(label)))

                        # Construct seeds from label boundary vertices
                        seeds = -1 * np.ones(len(points))
                        for ID, label_pair_list in enumerate(label_pair_lists):
                            label_pairs = [x for x in label_pair_list if label in x]
                            for label_pair in label_pairs:
                                indices_pair = [x for i,x in enumerate(indices_fold_pairs)
                                    if np.sort(fold_pairs[i]).tolist() == label_pair]
                                if len(indices_pair):

                                    # Do not include short boundary segments
                                    if min_boundary > 1:
                                        indices_pair2 = []
                                        seeds2 = segment(indices_pair, neighbor_lists)
                                        for seed2 in range(int(max(seeds2))+1):
                                            iseed2 = [i for i,x in enumerate(seeds2)
                                                      if x == seed2]
                                            if len(iseed2) >= min_boundary:
                                                indices_pair2.extend(iseed2)
                                            else:
                                                if len(iseed2) == 1:
                                                    print("    Remove assignment "
                                                          "of ID {0} from 1 vertex".
                                                          format(seed2))
                                                else:
                                                    print("    Remove assignment "
                                                          "of ID {0} from {1} vertices".
                                                          format(seed2, len(iseed2)))
                                        indices_pair = indices_pair2

                                    # Assign sulcus IDs to seeds
                                    seeds[indices_pair] = ID

                        # Identify vertices with the label
                        label_array = -1 * np.ones(len(points))
                        indices_label = [fold[i] for i,x in enumerate(fold_labels)
                                         if x == label]
                        if len(indices_label):
                            label_array[indices_label] = 1

                            # Propagate from seeds to vertices with label
                            #indices_seeds = []
                            #for seed in range(int(max(seeds))+1):
                            #    indices_seeds.append([i for i,x in enumerate(seeds)
                            #                          if x == seed])
                            #sulci2 = segment(indices_label, neighbor_lists,
                            #                 50, indices_seeds, False, True, labels)
                            sulci2 = propagate(points, faces, label_array,
                                               seeds, sulci, max_iters=10000,
                                               tol=0.001, sigma=5)
                            sulci[sulci2 > -1] = sulci2[sulci2 > -1]

    #---------------------------------------------------------------------------
    # Print out assigned sulci
    #---------------------------------------------------------------------------
    sulcus_numbers = [int(x) for x in np.unique(sulci) if x > -1]
    n_sulci = len(sulcus_numbers)
    print("Extracted {0} sulci from {1} folds ({2:.1f}s):".
          format(n_sulci, n_folds, time()-t0))
    if len(sulcus_names):
        for sulcus_number in sulcus_numbers:
            print("  {0}: {1}".format(sulcus_number, sulcus_names[sulcus_number]))
    else:
        print("  " + ", ".join([str(x) for x in sulcus_numbers]))

    #---------------------------------------------------------------------------
    # Print out unresolved sulci
    #---------------------------------------------------------------------------
    unresolved = [i for i in range(len(label_pair_lists))
                  if i not in sulcus_numbers]
    if len(unresolved) == 1:
        print("The following sulcus is unaccounted for:")
    else:
        print("The following {0} sulci are unaccounted for:".format(len(unresolved)))
    if len(sulcus_names):
        for sulcus_number in unresolved:
            print("  {0}: {1}".format(sulcus_number, sulcus_names[sulcus_number]))
    else:
        print("  " + ", ".join([str(x) for x in unresolved]))

    return sulci, n_sulci

#===============================================================================
# Example: extract_sulci()
#===============================================================================
if __name__ == "__main__":

    import os
    from time import time
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, write_vtk
    from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    from mindboggle.extract.extract_folds import extract_sulci
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    data_path = os.environ['MINDBOGGLE_DATA']

    # Load labels, folds, neighbor lists, and sulcus names and label pairs
    folds_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
                                        'features', 'lh.folds.vtk')
    labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
                              'labels', 'lh.labels.DKT25.manual.vtk')
    faces, lines, indices, points, npoints, folds, name = read_vtk(folds_file)
    faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file)
    neighbor_lists = find_neighbors(faces, npoints)
    sulcus_names_file = os.path.join(data_path, 'info', 'sulcus_names.txt')
    fid = open(sulcus_names_file, 'r')
    sulcus_names = fid.readlines()
    sulcus_names = [x.strip('\n') for x in sulcus_names]
    label_pair_lists = sulcus_boundaries()
    min_boundary = 1

    # Extract sulci
    sulci, n_sulci = extract_sulci(labels_file, folds, labels,
                                   neighbor_lists, label_pair_lists,
                                   min_boundary, sulcus_names)

    # Finally, write points, faces and sulci to a new vtk file
    #rewrite_vtk(labels_file, 'test_extract_sulci.vtk',
    #    [sulci.tolist()], ['sulci'], sulci.tolist())
    indices = [i for i,x in enumerate(sulci) if x > -1]
    write_vtk('test_extract_sulci.vtk', points, indices, lines,
              inside_faces(faces, indices), [sulci.tolist()], ['sulci'])
    os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')
