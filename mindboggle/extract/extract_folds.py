#!/usr/bin/python
"""
Functions to extract folds, sulci, or identify sulci from folds.

Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import sys
#import numpy as np
#from time import time
#from mindboggle.utils.io_vtk import load_scalars, load_points_faces,\
#    write_scalar_lists, rewrite_scalar_lists
#from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries,\
#    segment, fill_holes, compute_distance
#from mindboggle.info.sulcus_boundaries import sulcus_boundaries

#===============================================================================
# Extract folds
#===============================================================================
def extract_folds(depth_file, min_fold_size=1, do_fill_holes=False):
    """
    Use depth to extract folds from a triangular surface mesh.

    To extract an initial set of deep vertices from the surface mesh,
    we anticipate that there will be a rapidly decreasing distribution
    of low depth values (on the outer surface) with a long tail
    of higher depth values (in the folds), so we smooth the histogram's
    bin values (Gaussian), convolve to compute slopes,
    and find the depth value for the first bin with slope = 0.

    The resulting separately numbered folds may have holes
    resulting from shallower areas within a fold,
    so we call fill_holes(), which removes the largest region boundary,
    leaving smaller boundaries, presumably contours of holes within a region,
    and calls label_holes() to fill holes with surrounding region numbers.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    min_fold_size : int
        minimum fold size (number of vertices)
    do_fill_holes : Boolean
        fill holes?

    Returns
    -------
    folds : array of integers
        fold numbers for all vertices (default -1 for non-fold vertices)
    n_folds :  int
        number of folds

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalars, rewrite_scalar_lists
    >>> from mindboggle.utils.mesh_operations import inside_faces
    >>> from mindboggle.extract.extract_folds import extract_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>>
    >>> folds, n_folds = extract_folds(depth_file, 50, False)
    >>>
    >>> # Write results to vtk file and view with mayavi2:
    >>> folds = folds.tolist()
    >>> rewrite_scalar_lists(depth_file, 'test_extract_folds.vtk', [folds], ['folds'], folds)
    >>> #points, faces, depths, n_vertices = load_scalars(depth_file, True)
    >>> #indices = [i for i,x in enumerate(folds) if x > -1]
    >>> #write_scalar_lists('test_extract_folds.vtk', points, indices,
    >>> #    inside_faces(faces, indices), [folds], ['folds'])
    >>> os.system('mayavi2 -m Surface -d test_extract_folds.vtk &')

    """
    import numpy as np
    from time import time
    from scipy.ndimage.filters import gaussian_filter1d
    from mindboggle.utils.io_vtk import load_scalars
    from mindboggle.utils.mesh_operations import find_neighbors, segment, \
        fill_holes, fill_boundaries, propagate

    print("Extract folds in surface mesh")
    t0 = time()

    # Load depth and surface area values for all vertices
    points, faces, depths, n_vertices = load_scalars(depth_file, True)

    # Find neighbors for each vertex
    neighbor_lists = find_neighbors(faces, len(points))

    # Compute histogram of depth measures
    min_vertices = 10000
    if n_vertices > min_vertices:
        nbins = np.round(n_vertices / 100.0)
    else:
        error("  Expecting at least {0} vertices to create depth histogram".
        format(min_vertices))
    bins, bin_edges = np.histogram(depths, bins=nbins)
    """
    >>> # Plot histogram:
    >>> a,b,c = hist(depths, bins=nbins)
    """

    # Anticipating that there will be a rapidly decreasing distribution
    # of low depth values (on the outer surface) with a long tail of higher
    # depth values (in the folds), smooth the bin values (Gaussian), convolve
    # to compute slopes, and find the depth for the first bin with slope = 0.
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

    # Iteratively extract vertices from deep to less deep
    number_of_thresholds = 3
    thresholds = np.linspace(1, depth_threshold, num=number_of_thresholds+3)[3::]

    # Find the deepest vertices (depths greater than the highest depth threshold)
    for ithreshold in range(number_of_thresholds):
        indices_deep = [i for i,x in enumerate(depths) if x >= thresholds[ithreshold]]
        if len(indices_deep):
            break
    if len(indices_deep):

        # Segment initial set of folds
        print("  Segment vertices deeper than {0:.2f}".format(thresholds[ithreshold]))
        t1 = time()
        folds = segment(indices_deep, neighbor_lists)
        # Slightly slower alternative -- fill boundaries:
        #regions = -1 * np.ones(len(points))
        #regions[indices_deep] = 1
        #folds = fill_boundaries(regions, neighbor_lists)
        print('    ...Segmented deepest vertices ({0:.2f} seconds)'.format(time() - t1))
        """
        >>> # Display resulting initial folds:
        >>> rewrite_scalar_lists(depth_file, 'test_folds1.vtk', [folds.tolist()], 'folds', folds)
        >>> os.system('mayavi2 -m Surface -d test_folds1.vtk &')
        """

        # Expand folds iteratively
        print("  Grow folds by including shallower vertices")
        for threshold in thresholds[ithreshold+1::]:

            indices_deep = [i for i,x in enumerate(depths) if x >= threshold]
            unique_folds = [x for x in np.unique(folds) if x > -1]
            fold_lists = [[] for x in unique_folds]
            for ifold, nfold in enumerate(unique_folds):
                fold_lists[ifold] = [i for i,x in enumerate(folds) if x == nfold]
            folds2 = segment(indices_deep, neighbor_lists, 1,
                             fold_lists, keep_seeding=True)
            folds[folds2 > -1] = folds2[folds2 > -1]
            #folds = propagate(points, faces, deep_vertices, folds, folds,
            #                  max_iters=10000, tol=0.001, sigma=5)

        print('    ...Segmented folds ({0:.2f} seconds)'.format(time() - t1))
        n_folds = len([x for x in list(set(folds)) if x != -1])

        # If there are any folds, find and fill holes
        if n_folds > 0 and do_fill_holes:
            folds = fill_holes(folds, neighbor_lists)

        # Remove small folds
        if min_fold_size > 1:
            print('    Remove folds smaller than {0}'.format(min_fold_size))
            for nfold in np.unique(folds):
                indices_fold = np.where(folds == nfold)[0]
                if len(indices_fold) < min_fold_size:
                    folds[indices_fold] = -1
            n_folds = len(np.unique(folds))

        print('  ...Extracted {0} folds ({0:.2f} seconds)'.
              format(n_folds, time() - t0))
    else:
        print('  No deep vertices')

    # Return folds, number of folds
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

        A ``fold`` is a group of connected, deep vertices.

        A ``label list`` is the list of labels used to define a single sulcus.

        A ``label pair list`` contains pairs of labels, where each pair
        defines a boundary between two labeled regions.
        No two label pair lists share a label pair.

        A ``sulcus ID`` uniquely identifies a sulcus, as index to a sulcus label pair list.

    Algorithm ::

      For each fold (vertices with the same fold number):

        **Case 0**: NO MATCH -- fold has no label

        **Case 1**: NO MATCH -- fold has only one label

          If the fold has only one label, remove the fold by assigning
          -1 to all vertices in this fold.

        **Case 2**: matching label set in sulcus

          If the set of labels in the fold is the same as in one of the
          protocol's sulci, assign the index to the corresponding
          label list to all vertices in the fold.

        **Case 3**: NO MATCH -- fold has no sulcus label pair

        **Case 4**: fold labels in one sulcus

          If the set of labels in the fold is a subset of labels in
          only one of the protocol's sulcus label lists, assign the index
          for that list to all vertices in the fold.
          (Ex: the labels in the fold are [1,2,3] and there is only one
           sulcus label list containing 1, 2, and 3: [1,2,3,4])

        **Case 5**: ambiguous -- fold labels in more than one sulcus
                              -- fold labels not contained by a sulcus

        **Case 6**: vertex labels in only one of the fold's sulcus label pairs

          Find vertices with labels that are in only one of the fold's
          label boundary pairs. Assign the vertices the sulcus with the
          label pair if they are connected to the label boundary for that pair,
          via label propagation or seed growing.

        **Case 7**: remaining vertices connected to sulcus label boundaries

          If there are remaining vertices, segment into sets of vertices
          connected to label boundary seeds (remaining label boundary vertices),
          and assign a sulcus ID to each segment.

    Parameters
    ----------
    surface_vtk : string
        file name for surface mesh vtk (from which to extract points and faces)
    labels : list of integers
        labels for all vertices
    fold_lists : list of lists of integers
        each list contains indices to vertices of a fold
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
    >>> from mindboggle.utils.io_vtk import load_points_faces
    >>> from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    >>> from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    >>> from mindboggle.extract.extract_folds import extract_sulci
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>>
    >>> # Load labels, folds, neighbor lists, and sulcus names and label pairs
    >>> folds_file = os.path.join(data_path, 'results', 'features',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'folds.vtk')
    >>> labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> points, faces, folds, n_vertices = load_scalars(folds_file, False)
    >>> points, faces, labels, n_vertices = load_scalars(labels_file, False)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> sulcus_names_file = os.path.join(data_path, 'info', 'sulcus_names.txt')
    >>> fid = open(sulcus_names_file, 'r')
    >>> sulcus_names = fid.readlines()
    >>> sulcus_names = [x.strip('\n') for x in sulcus_names]
    >>> label_pair_lists = sulcus_boundaries()
    >>> min_boundary = 10
    >>> vtk_file = 'test_extract_sulci.vtk'
    >>>
    >>> # Extract sulci
    >>> sulci, n_sulci = extract_sulci(labels_file, folds, labels,
    >>>                                neighbor_lists, label_pair_lists,
    >>>                                min_boundary, sulcus_names)
    >>>
    >>> # Finally, write points, faces and sulci to a new vtk file
    >>> #rewrite_scalar_lists(labels_file, vtk_file,
    >>> #                     [sulci.tolist()], ['sulci'], sulci.tolist())
    >>> indices = [i for i,x in enumerate(sulci) if x > -1]
    >>> write_scalar_lists('test_extract_sulci.vtk', points, indices,
    >>>    inside_faces(faces, indices), [sulci.tolist()], ['sulci'])
    >>> os.system('mayavi2 -m Surface -d ' + vtk_file + ' &')

    """
    from time import time
    import numpy as np
    from mindboggle.utils.io_vtk import load_points_faces
    from mindboggle.utils.mesh_operations import detect_boundaries, propagate, segment
    from mindboggle.label.label_functions import find_superset_subset_lists

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
    points, faces = load_points_faces(surface_vtk)

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    fold_numbers = [int(x) for x in np.unique(folds) if x > -1]
    n_folds = len(fold_numbers)
    print("Extract sulci from {0} folds...".format(n_folds))
    t0 = time()
#    for n_fold in fold_numbers:
    for n_fold in [1]:

        fold = np.where(folds == n_fold)[0]
        len_fold = len(fold)
        ambiguous_case = False

        # List the labels in this fold (greater than zero)
        fold_labels = [labels[x] for x in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels) if x > 0]

        #-----------------------------------------------------------------------
        # Case 0: NO MATCH -- fold has no label
        #-----------------------------------------------------------------------
        if not len(unique_fold_labels):
            print("  Fold {0} of {1} ({2} vertices): "
                  "NO MATCH -- fold has no label".
                  format(n_fold + 1, n_folds, len_fold))
            # Ignore: sulci already initialized with -1 values

        #-----------------------------------------------------------------------
        # Case 1: NO MATCH -- fold has only one label
        #-----------------------------------------------------------------------
        elif len(unique_fold_labels) == 1:
            print("  Fold {0} of {1} ({2} vertices): "
                  "NO MATCH -- fold has only one label ({3})".
                  format(n_fold + 1, n_folds, len_fold, unique_fold_labels[0]))
            # Ignore: sulci already initialized with -1 values

        #-----------------------------------------------------------------------
        # Case 2: matching label set in sulcus
        #-----------------------------------------------------------------------
        elif unique_fold_labels in label_lists:
            if len(sulcus_names):
                print("  Fold {0} of {1} ({2} vertices): matching label set "
                      "in sulcus: {3} ({4})".format(n_fold + 1, n_folds, len_fold,
                           sulcus_names[label_lists.index(unique_fold_labels)],
                           ', '.join([str(x) for x in unique_fold_labels])))
            else:
                print("  Fold {0} of {1} ({2} vertices): matching label set "
                      "in sulcus ({3})".format(n_fold + 1, n_folds, len_fold,
                           ', '.join([str(x) for x in unique_fold_labels])))
            # Assign ID of the matching sulcus to all fold vertices
            sulci[fold] = label_lists.index(unique_fold_labels)

        # Cases 3-5: at least one fold label but no perfect match
        else:

            # Find all label boundary pairs within the fold
            indices_fold_pairs, fold_pairs, unique_fold_pairs = detect_boundaries(
                fold, labels, neighbor_lists)

            # Find fold label pairs in the protocol (pairs are already sorted)
            fold_pairs_in_protocol = [x for x in unique_fold_pairs
                                      if x in protocol_pairs]

            print("  Fold {0} labels: {1}".format(n_fold + 1,
                  ', '.join([str(x) for x in unique_fold_labels])))
            print("  Fold {0} label pairs in protocol: {1}".format(n_fold+1,
                  ', '.join([str(x) for x in fold_pairs_in_protocol])))

            #-------------------------------------------------------------------
            # Case 3: NO MATCH -- fold has no sulcus label pair
            #-------------------------------------------------------------------
            if not len(fold_pairs_in_protocol):
                print("  Fold {0} of {1} ({2} vertices): "
                      "NO MATCH -- fold has no sulcus label pair".
                      format(n_fold + 1, n_folds, len_fold))

            # Cases 4-5: labels in common between fold and sulcus/sulci
            else:
                # Find overlap of sulcus labels and fold labels
                superset_indices, subset_indices = find_superset_subset_lists(
                    unique_fold_labels, label_lists)

                #---------------------------------------------------------------
                # Cases 4: fold labels contained by one sulcus
                #---------------------------------------------------------------
                if len(superset_indices) == 1:
                    if len(sulcus_names):
                        print("  Fold {0} of {1} ({2} vertices): "
                              "fold labels in one sulcus: {3} ({4})".
                        format(n_fold + 1, n_folds, len_fold,
                               sulcus_names[superset_indices[0]],
                               ', '.join([str(x)
                                    for x in label_lists[superset_indices[0]]])))
                    else:
                        print("  Fold {0} of {1} ({2} vertices): "
                              "fold labels in one sulcus: ({3})".
                            format(n_fold + 1, n_folds, len_fold,
                            ', '.join([str(x)
                                 for x in label_lists[superset_indices[0]]])))
                    # Assign ID of matching sulcus to all fold vertices
                    sulci[fold] = superset_indices[0]

                #---------------------------------------------------------------
                # Case 5: ambiguous  -- fold labels in more than one sulcus
                #                    -- fold labels not contained by a sulcus
                #---------------------------------------------------------------
                else:
                    print("  Fold {0} of {1} ({2} vertices): ambiguous -- "
                          "fold labels contained by multiple or by no sulci".
                          format(n_fold + 1, n_folds, len_fold))
                    ambiguous_case = True

        #-----------------------------------------------------------------------
        # Ambiguous case
        #-----------------------------------------------------------------------
        if ambiguous_case:

            # Labels in the protocol (includes repeats across label pairs)
            labels_in_pairs = [x for lst in fold_pairs_in_protocol for x in lst]

            # Labels that appear in one or in more than one sulcus label boundary
            unique_labels = []
            nonunique_labels = []
            for label in np.unique(labels_in_pairs):
                if len([x for x in labels_in_pairs if x == label]) == 1:
                    unique_labels.append(label)
                else:
                    nonunique_labels.append(label)

            #-------------------------------------------------------------------
            # Case 6: vertices whose labels are in only one sulcus label pair
            #-------------------------------------------------------------------
            # Find vertices with a label that is in only one of the fold's
            # label pairs (the other label in the pair can exist
            # in other pairs). Assign the vertices the sulcus with the label
            # pair if they are connected to the label boundary for that pair.
            #-------------------------------------------------------------------
            if len(unique_labels):
                for pair in fold_pairs_in_protocol:

                    # If one or both labels in label pair is/are unique
                    unique_labels_in_pair = [x for x in pair if x in unique_labels]
                    n_unique = len(unique_labels_in_pair)
                    if len(unique_labels_in_pair):

                        ID = [i for i,x in enumerate(label_pair_lists) if pair in x][0]

                        # Construct seeds from label boundary vertices
                        indices_pair = [x for i,x in enumerate(indices_fold_pairs)
                                        if list(set(fold_pairs[i])) == pair]

                        # Identify vertices with unique label(s) in pair
                        indices_unique_labels = [fold[i]
                                                 for i,x in enumerate(fold_labels)
                                                 if x in unique_labels_in_pair]

                        # Propagate from seeds to labels in label pair
                        sulci2 = segment(indices_unique_labels, neighbor_lists,
                                         min_region_size=1,
                                         seed_lists=[indices_pair],
                                         keep_seeding=False,
                                         spread_within_labels=True, labels=labels)
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

            #-------------------------------------------------------------------
            # Case 7: vertex labels shared by multiple label pairs
            #-------------------------------------------------------------------
            # Propagate labels from label boundaries to vertices with labels
            # that are shared by multiple label pairs in the fold.
            #-------------------------------------------------------------------
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
                    label_array[indices_label] = 1

                    # Propagate from seeds to vertices with label
                    #indices_seeds = []
                    #for seed in range(int(max(seeds))+1):
                    #    indices_seeds.append([i for i,x in enumerate(seeds)
                    #                          if x == seed])
                    #sulci2 = segment(indices_label, neighbor_lists,
                    #                 50, indices_seeds, False, True, labels)
                    sulci2 = propagate(points, faces, label_array, seeds, sulci,
                                       max_iters=10000, tol=0.001, sigma=5)
                    sulci[sulci2 > -1] = sulci2[sulci2 > -1]

    # Print out assigned sulci
    sulcus_numbers = [int(x) for x in np.unique(sulci) if x > -1]
    n_sulci = len(sulcus_numbers)
    print("Extracted {0} sulci from {1} folds ({2:.1f}s):".
          format(n_sulci, n_folds, time()-t0))
    if len(sulcus_names):
        for sulcus_number in sulcus_numbers:
            print("  {0}: {1}".format(sulcus_number, sulcus_names[sulcus_number]))
    else:
        print("  " + ", ".join([str(x) for x in sulcus_numbers]))

    # Print out unresolved sulci
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


# Extract_sulci()
if __name__ == "__main__":

    import os
    from time import time
    import numpy as np
    from mindboggle.utils.io_vtk import load_points_faces
    from mindboggle.utils.io_vtk import load_scalars, write_scalar_lists
    from mindboggle.utils.mesh_operations import find_neighbors, inside_faces
    from mindboggle.extract.extract_folds import extract_sulci
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    data_path = os.environ['MINDBOGGLE_DATA']

    # Load labels, folds, neighbor lists, and sulcus names and label pairs
    folds_file = os.path.join(data_path, 'results', 'features',
                '_hemi_lh_subject_MMRR-21-1', 'folds.vtk')
    labels_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
                 'label', 'lh.labels.DKT25.manual.vtk')
    points, faces, folds, n_vertices = load_scalars(folds_file, True)
    points, faces, labels, n_vertices = load_scalars(labels_file, True)
    neighbor_lists = find_neighbors(faces, len(points))
    sulcus_names_file = os.path.join(data_path, 'info', 'sulcus_names.txt')
    fid = open(sulcus_names_file, 'r')
    sulcus_names = fid.readlines()
    sulcus_names = [x.strip('\n') for x in sulcus_names]
    label_pair_lists = sulcus_boundaries()

    # Extract sulci
    sulci, n_sulci = extract_sulci(labels_file, folds, labels, neighbor_lists,
                                   label_pair_lists, min_boundary=10,
                                   sulcus_names=sulcus_names)

    # Finally, write points, faces and sulci to a new vtk file
    #rewrite_scalar_lists(labels_file, 'test_extract_sulci.vtk',
    #    [sulci.tolist()], ['sulci'], sulci.tolist())
    indices = [i for i,x in enumerate(sulci) if x > -1]
    write_scalar_lists('test_extract_sulci.vtk', points, indices,
        inside_faces(faces, indices), [sulci.tolist()], ['sulci'])
    os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')
