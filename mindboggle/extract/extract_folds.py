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
#from mindboggle.utils.io_vtk import load_scalar, write_scalars, rewrite_scalars
#from mindboggle.utils.mesh_operations import find_neighbors, detect_boundaries,\
#    segment, fill_holes, compute_distance
#from mindboggle.info.sulcus_boundaries import sulcus_boundaries

#===============================================================================
# Extract folds
#===============================================================================
def extract_all_folds(depth_file, area_file, fraction_folds):
    """
    Use depth to extract all folds from a triangular surface mesh,
    by fraction of surface area.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    area_file : str
        surface mesh file in VTK format with faces and surface area scalar values
    fraction_folds : float
        fraction of surface mesh considered folds

    Returns
    -------
    folds : array of integers
        an integer for every mesh vertex: 1 for fold, -1 for non-fold

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, write_scalars
    >>> from mindboggle.extract.extract_folds import extract_all_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> folds, min_depth, n_deep_vertices = extract_all_folds(depth_file,
    >>>     area_file, 0.5)
    >>> # Write results to vtk file and view with mayavi2:
    >>> from mindboggle.utils.io_vtk import rewrite_scalars
    >>> rewrite_scalars(depth_file, 'test_extract_all_folds.vtk', folds, folds)
    >>> os.system('mayavi2 -m Surface -d test_extract_all_folds.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.utils.io_vtk import load_scalar

    print("Extract the deepest surface mesh vertices ({0} of surface area)...".
          format(fraction_folds))
    t0 = time()

    # Load depth and surface area values from VTK files
    points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=True)
    points, faces, areas, n_vertices = load_scalar(area_file, return_arrays=True)

    sort_indices = np.argsort(depths)
    sort_indices = sort_indices[::-1]

    total_area = np.sum(areas)
    fraction_area = fraction_folds * total_area
    sum_area = 0
    folds = -1 * np.ones(len(areas))
    min_depth = 1
    for index in sort_indices:
        folds[index] = 1
        sum_area += areas[index]
        if sum_area >= fraction_area:
            min_depth = depths[index]
            break

    n_deep_vertices = len([x for x in folds if x == 1])

    print('  ...Extracted {0} vertices deeper than {1:.2f} ({2:.2f} seconds)'.
          format(n_deep_vertices, min_depth, time() - t0))

    return folds, min_depth, n_deep_vertices

#===============================================================================
# Extract individual folds
#===============================================================================
def extract_folds(depth_file, area_file, neighbor_lists, fraction_folds,
                  min_fold_size, do_fill_holes=False):
    """
    Use depth to extract folds from a triangular surface mesh and fill holes
    resulting from shallower areas within a fold.

    Parameters
    ----------
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    area_file : str
        surface mesh file in VTK format with faces and surface area scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    fraction_folds : float
        fraction of surface mesh considered folds
    min_fold_size : int
        minimum fold size (number of vertices)
    do_fill_holes : Boolean
        segment and fill holes?

    Returns
    -------
    fold_IDs : array of integers
        fold IDs for all vertices, with -1s for non-fundus vertices
    n_folds :  int
        number of folds

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.extract.extract_folds import extract_folds
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> points, faces, depths, n_vertices = load_scalar(depth_file, return_arrays=0)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> fold_IDs, n_folds = extract_folds(depth_file, area_file, neighbor_lists, 0.5, 50, False)
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalars(depth_file, 'test_extract_folds.vtk', fold_IDs, fold_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_folds.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.utils.mesh_operations import segment, fill_holes
    from mindboggle.utils.io_vtk import load_scalar
    from mindboggle.extract.extract_folds import extract_all_folds

    print("Extract folds from surface mesh...")
    t0 = time()

    # Compute the minimum depth threshold for defining folds
    folds, min_depth, n_deep_vertices = extract_all_folds(depth_file,
        area_file, fraction_folds)

    # Segment folds of a surface mesh
    print("  Segment surface mesh into separate folds deeper than {0:.2f}...".
          format(min_depth))
    t1 = time()
    vertices_to_segment = [i for i,x in enumerate(folds) if x == 1]
    fold_IDs = segment(vertices_to_segment, neighbor_lists,
        seed_lists=[], min_region_size=min_fold_size,
        spread_same_labels=False, labels=[], label_pair_lists=[])
    print('    ...Folds segmented ({0:.2f} seconds)'.format(time() - t1))
    n_folds = len([x for x in list(set(fold_IDs)) if x != -1])

    # If there are any folds
    if n_folds > 0 and do_fill_holes:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep)
        t2 = time()
        vertices_to_segment = [i for i,x in enumerate(fold_IDs) if x==-1]

        # Segment holes in the folds
        print('  Segment holes in the folds...')
        holes = segment(vertices_to_segment, neighbor_lists,
                seed_lists=[], min_region_size=1,
                spread_same_labels=False, labels=[], label_pair_lists=[])
        n_holes = len([x for x in list(set(holes))])

        # If there are any holes
        if n_holes > 0:

            # Ignore the largest hole (the background) and renumber holes
            max_hole_size = 0
            max_hole_index = 0
            for ihole in range(n_holes):
                I = np.where(holes == ihole)
                if len(I) > max_hole_size:
                    max_hole_size = len(I)
                    max_hole_index = ihole
            holes[holes == max_hole_index] = -1
            if max_hole_index < n_holes:
                holes[holes > max_hole_index] -= 1
            n_holes -= 1
            print('    ...{0} holes segmented ({1:.2f} seconds)'.
                  format(n_holes, time() - t2))

            # Fill holes
            t3 = time()
            fold_IDs = fill_holes(fold_IDs, holes, n_holes, neighbor_lists)
            print('  Filled holes ({0:.2f} seconds)'.format(time() - t3))

    print('  ...Extracted folds greater than {0:.2f} depth in {1:.2f} seconds'.
          format(min_depth, time() - t0))

    # Return folds, number of folds
    return fold_IDs, n_folds

#===============================================================================
# Identify sulci from folds
#===============================================================================
def identify_sulci_from_folds(labels, folds, neighbor_lists, sulcus_names,
                              label_pair_lists):
    """
    Identify sulcus folds in a brain surface according to a labeling protocol
    that includes a list of label pairs defining each sulcus.

    Parameters
    ----------
    labels : list of integers
        labels for all vertices
    folds : list of lists of integers
        each list contains indices to vertices of a fold
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices for each vertex
    sulcus_names : list of strings
        names of sulci
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus

    Returns
    -------
    sulcus_IDs : array of integers
        sulcus IDs for all vertices, with -1s for non-sulcus vertices

    Definitions
    -----------

    A fold is a group of connected vertices deeper than a depth threshold.

    A ''label list'' is the list of labels used to define a single sulcus.

    A ''label pair list'' contains pairs of labels, where each pair
    defines a boundary between two labeled regions.
    No two label pair lists share a label pair.

    A ''sulcus ID'' uniquely identifies a sulcus.
    It is the index to a sulcus label list (or sulcus label pair list).

    Algorithm
    ---------

    First, remove all fold vertices with labels that are not in the label list
    by assigning -1 to them.

    For each fold (vertices with the same non-zero label):

        Case 0: NO MATCH -- fold has no label

        Case 1: NO MATCH -- fold has only one label

          If the fold has only one label, remove the fold by assigning a
          sulcus ID of -1 to all vertices in this fold.

        Case 2: matching label set in sulcus

          If the set of labels in the fold is the same as in one of the
          protocol's sulci, assign the index to the corresponding sulcus
          label list to all vertices in the fold.

        Case 3: NO MATCH -- fold has no sulcus label pair

        Case 4a: fold labels in one sulcus

          If the set of labels in the fold is a subset of labels in
          only one of the protocol's sulcus label lists, assign the index
          for that list to all vertices in the fold.
          (Ex: the labels in the fold are [1,2,3] and there is only one
           sulcus label list containing 1, 2, and 3: [1,2,3,4])

        Case 4b: labels in sulci but label pair in one sulcus

          If the set of labels in the fold is a subset of labels in more than
          one of the protocol's sulcus label lists, find corresponding sulcus
          label pair lists that share a label pair with the fold label pairs.

        Case 4c: ambiguous -- fold label pairs in multiple sulci

        Case 5: ambiguous -- fold labels not contained by a sulcus

        For ambiguous cases above:

          Find label boundary pairs in the fold whose labels
          are shared by any other label pairs in the fold,
          and store the sulcus IDs for these pairs.

          Assign a sulcus ID to fold vertices that have unique
          label pair labels (unique to ensure a label is not shared
          with another label pair); store unassigned vertices.

          If there are remaining vertices with duplicate label pair labels,
          construct seed lists of remaining label boundary vertices,
          segment into sets of vertices connected to label boundary seeds,
          and assign a sulcus ID to each segment.

    import sys
    from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    from mindboggle.utils.io_file import read_columns
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    from mindboggle.utils.mesh_operations import find_neighbors

    Example
    -------
    >>> # Arguments
    >>> labels_file = sys.argv[1]
    >>> folds_file = sys.argv[2]
    >>> sulcus_names_file = sys.argv[3]
    >>> vtk_file = sys.argv[4]
    >>> # Load labels and folds (the second surface has to be inflated).
    >>> points, faces, labels, n_vertices = load_scalar(labels_file, return_arrays=0)
    >>> points, faces, fold_IDs, n_vertices = load_scalar(folds_file, return_arrays=0)
    >>> fid = open(sulcus_names_file, 'r')
    >>> sulcus_names = fid.readlines()
    >>> sulcus_names = [x.strip('\n') for x in sulcus_names]
    >>> # Calculate neighbor lists for all vertices
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> # Prepare list of all unique sorted label pairs in the labeling protocol
    >>> label_pair_lists = sulcus_boundaries()
    >>> # Create a list of lists of folds
    >>> unique_fold_IDs = set(fold_IDs)
    >>> folds = []
    >>> for id in unique_fold_IDs:
    >>>     if id > 0:
    >>>         fold = [i for i,x in enumerate(fold_IDs) if x == id]
    >>>         folds.append(fold)
    >>> # Identify sulci from folds
    >>> sulcus_IDs = identify_sulci_from_folds(labels, folds, neighbor_lists,
    >>>                       sulcus_names, label_pair_lists)
    >>> # Finally, write points, faces and sulcus_IDs to a new vtk file
    >>> rewrite_scalars(labels_file, vtk_file, sulcus_IDs, filter_scalars=sulcus_IDs)

    """
    import numpy as np
    from mindboggle.utils.mesh_operations import segment, detect_boundaries, compute_distance

    verbose = 0

    #---------------------------------------------------------------------------
    # Nested function definitions
    #---------------------------------------------------------------------------
    def find_superset_subset_lists(labels, label_lists):
        """
        Find *label_lists* that are supersets or subsets of *labels*.

        Parameters
        ----------
        labels : list of integers
            label numbers
        label_lists : list of lists of integers
            each list contains label numbers

        Returns
        -------
        superset_indices : list of integers
            indices to label_lists that are a superset of labels
        subset_indices : list of integers
            indices to label_lists that are a subset of labels

        Example
        -------
        >>> find_superset_subset_lists([1,2],[[1,2],[3,4]])
        >>> [0]
        >>> find_superset_subset_lists([1,2],[[1,2,3],[1,2,5]])
        >>> [0, 1]
        >>> find_superset_subset_lists([1,2],[[2,3],[1,2,5]])
        >>> [1]

        """

        labels = set(labels)
        superset_indices = []
        subset_indices = []
        for Id, label_list in enumerate(label_lists):
            if labels.issubset(set(label_list)):
                superset_indices.append(Id)
            if set(label_list).issubset(labels):
                subset_indices.append(Id)

        return superset_indices, subset_indices

    #---------------------------------------------------------------------------
    # Prepare data structures
    #---------------------------------------------------------------------------
    # Array of sulcus IDs for fold vertices, initialized as -1.
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in the label list, or vertices having only one label,
    # their sulcus IDs will remain -1.
    sulcus_IDs = -1 * np.ones(len(neighbor_lists))

    # Prepare list of sulcus label lists (labels within a sulcus)
    label_lists = []
    for row in label_pair_lists:
        label_lists.append(list(np.unique(np.asarray(row))))

    # Prepare list of all unique sorted label pairs in the labeling protocol
    protocol_pairs = []
    [protocol_pairs.append(list(set(x)))
     for lst in label_pair_lists
     for x in lst
     if list(set(x)) not in protocol_pairs]

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    for ifold, fold in enumerate(folds):
        special_case = False

        # List the labels in this fold
        fold_labels = [labels[vertex] for vertex in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels)]

        # Case 0: NO MATCH -- fold has no label
        if not len(unique_fold_labels):
            print("  Fold {0}: NO MATCH -- fold has no label")
            # Ignore: sulcus_IDs already initialized with -1 values
            continue

        # Case 1: NO MATCH -- fold has only one label
        elif len(unique_fold_labels) == 1:
            print("  Fold {0}: NO MATCH -- fold has only one label ({1})".
                  format(ifold, unique_fold_labels[0]))
            # Ignore: sulcus_IDs already initialized with -1 values
            continue

        # Case 2: matching label set in sulcus
        elif unique_fold_labels in label_lists:
            print("  Fold {0}: matching label set in sulcus: {1} ({2})".
                format(ifold,
                sulcus_names[label_lists.index(unique_fold_labels)],
                  ', '.join([str(x) for x in unique_fold_labels])))
            index = label_lists.index(unique_fold_labels)
            # Assign ID of the matching sulcus to all fold vertices
            sulcus_IDs[fold] = index

        # Cases 3-5: at least one fold label but no perfect match
        else:
            # Find all label boundary pairs within the fold
            indices_fold_pairs, fold_pairs, unique_fold_pairs = detect_boundaries(
                fold, labels, neighbor_lists)

            # Find fold label pairs in the protocol (pairs are already sorted)
            fold_pairs_in_protocol = [x for x in unique_fold_pairs
                                      if x in protocol_pairs]
            if verbose:
                print("  Fold labels: {0}".
                      format(', '.join([str(x) for x in unique_fold_labels])))
                #print("  Fold label pairs: {0}".
                #      format(', '.join([str(x) for x in unique_fold_pairs])))
                print("  Fold label pairs in protocol: {0}".
                      format(', '.join([str(x) for x in fold_pairs_in_protocol])))

            # Case 3: NO MATCH -- fold has no sulcus label pair
            if not len(fold_pairs_in_protocol):
                print("  Fold {0}: NO MATCH -- fold has no sulcus label pair".
                      format(ifold))

            # Cases 4-5: label pair(s) in common between sulci and fold
            else:

                # Find overlap of sulcus labels and fold labels
                superset_indices, subset_indices = find_superset_subset_lists(
                    unique_fold_labels, label_lists)

                # Cases 4a-c: fold labels contained by sulcus label list(s)
                if len(superset_indices):

                    # Case 4a: fold labels in one sulcus
                    if len(superset_indices) == 1:
                        print("  Fold {0}: fold labels in one sulcus: {1} ({2})".
                            format(ifold, sulcus_names[superset_indices[0]],
                            ', '.join([str(x)
                            for x in label_lists[superset_indices[0]]])))
                        # Assign ID of matching sulcus to all fold vertices
                        sulcus_IDs[fold] = superset_indices[0]

                    # Cases 4b-c: Fold labels in more than one sulcus
                    else:
                        # Find sulci that contain all of the fold's labels
                        # and share one or more label pairs with the fold
                        label_pair_array = np.array(label_pair_lists)
                        IDs = [superset_indices[i] for i, sublst
                               in enumerate(label_pair_array[superset_indices])
                               for x in sublst if x in fold_pairs_in_protocol]
                        IDs = list(set(IDs))

                        # Case 4b: labels in sulci but label pair in one sulcus
                        if len(IDs) == 1:
                            print("  Fold {0}: labels in sulci "
                                  "but label pair in one sulcus".format(ifold))
                            # Assign ID of matching sulcus to all fold vertices
                            sulcus_IDs[fold] = IDs

                        # Case 4c: ambiguous -- fold label pairs in multiple sulci
                        else:
                            print("  Fold {0}: ambiguous -- "
                                  "fold label pairs in multiple sulci".
                                  format(ifold))
                            special_case = True

                # Case 5: ambiguous -- fold labels not contained by a sulcus
                else:
                    print("  Fold {0}: ambiguous -- "
                          "fold labels not contained by a sulcus".format(ifold))
                    special_case = True

        # Special cases
        if special_case:

            # Find label boundary pairs in the fold whose labels
            # are and are not shared by any other label pairs
            # in the fold, and store the sulcus IDs for these pairs
            unique_pairs = []
            IDs_unique_pairs = []
            remainder_pairs = []
            IDs_remainder_pairs = []
            for pair in fold_pairs_in_protocol:
                labels_in_pairs = [x for sublst in fold_pairs_in_protocol
                                   for x in sublst]
                if len([x for x in labels_in_pairs if x in pair]) == 2:
                    unique_pairs.append(pair)
                    IDs_unique_pairs.extend(
                        [i for i,x in enumerate(label_pair_lists)
                         if np.sort(pair).tolist() in x])
                else:
                    remainder_pairs.append(pair)
                    IDs_remainder_pairs.extend([i
                        for i,x in enumerate(label_pair_lists)
                        if np.sort(pair).tolist() in x])

            # Assign a sulcus ID to fold vertices that have unique
            # label pair labels (unique to ensure a label is not shared
            # with another label pair); store unassigned vertices
            if len(unique_pairs):
                print("           Assign sulcus IDs to vertices "
                      "with labels in only one fold pair")
                if verbose:
                    for index, ID_unique_pair \
                      in enumerate(IDs_unique_pairs):
                        print("           {0} ({1})".format(
                            unique_pairs[index], sulcus_names[ID_unique_pair]))
                unassigned = []
                for vertex in fold:
                    index = [i for i,x in enumerate(unique_pairs)
                             if labels[vertex] in x]
                    if len(index):
                        sulcus_IDs[vertex] = IDs_unique_pairs[index[0]]
                    else:
                        unassigned.append(vertex)
            else:
                unassigned = fold[:]

            # If there are remaining vertices with duplicate label pair labels
            if len(unassigned) and len(remainder_pairs):

                # Construct seed lists of remaining label boundary vertices
                seed_lists = []
                for remainder_pair in remainder_pairs:
                    seed_lists.append([x for i,x in enumerate(indices_fold_pairs)
                        if list(set(fold_pairs[i])) == list(set(remainder_pair))])

                # Segment into sets of vertices connected to label boundary seeds
                print("           Segment into separate label-pair regions")
                subfolds = segment(unassigned, neighbor_lists,
                    seed_lists, min_region_size=50,
                    spread_same_labels=False, labels=[], label_pair_lists=[])
                #n_subfolds = len([x for x in list(set(subfolds)) if x != -1])

                # Assign a sulcus ID to each segment
                print("           Assign sulcus IDs to segmented vertices")
                for ifold in np.unique(subfolds).tolist():
                    if ifold > -1:
                        subfold = [i for i,x in enumerate(subfolds) if x == ifold]
                        sulcus_IDs[subfold] = IDs_remainder_pairs[ifold]

    return sulcus_IDs

#===============================================================================
# Extract sulci
#===============================================================================
def extract_sulci(label_pair_lists, labels, depth_file, area_file, neighbor_lists,
                  fraction_folds, min_sulcus_size, do_fill_holes=False):
    """
    Use depth and a sulcus labeling protocol to extract sulcus folds
    from an anatomically labeled triangular surface mesh and fill holes
    resulting from shallower areas within a fold.

    Parameters
    ----------
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    labels : array of integers
        label IDs for all vertices, with -1s for unlabeled vertices
    depth_file : str
        surface mesh file in VTK format with faces and depth scalar values
    area_file : str
        surface mesh file in VTK format with faces and surface area scalar values
    neighbor_lists : list of lists of integers
        each list contains indices to neighboring vertices
    fraction_folds : float
        fraction of surface mesh considered folds
    min_sulcus_size : int
        minimum sulcus size (number of vertices)
    do_fill_holes : Boolean
        segment and fill holes?

    Returns
    -------
    sulcus_IDs : array of integers
        sulcus IDs for all vertices, with -1s for non-sulcus vertices
    n_sulci :  int
        number of sulcus folds

    Definitions
    -----------

    A fold is a group of connected vertices deeper than a depth threshold.

    A ''label list'' is the list of labels used to define a single sulcus.

    A ''label pair list'' contains pairs of labels, where each pair
    defines a boundary between two labeled regions.
    No two label pair lists share a label pair.

    A ''sulcus ID'' uniquely identifies a sulcus.
    It is the index to a sulcus label list (or sulcus label pair list).

    Algorithm
    ---------

    First, extract folds by a depth threshold, and remove all non-fold
    vertices by assigning -1 to them.

    Construct seed lists of label boundary vertices, segment into sets
    of vertices connected to label boundary seeds that have the same labels,
    and assign a sulcus ID to each segment.

    Example
    -------
    >>> import os
    >>> from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    >>> from mindboggle.utils.mesh_operations import find_neighbors
    >>> from mindboggle.extract.extract_folds import extract_sulci
    >>> from mindboggle.info.sulcus_boundaries import sulcus_boundaries
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>              'label', 'lh.labels.DKT25.manual.vtk')
    >>> depth_file = os.path.join(data_path, 'measures',
    >>>              '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    >>> area_file = os.path.join(data_path, 'measures',
    >>>             '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')
    >>> points, faces, labels, n_vertices = load_scalar(label_file, True)
    >>> neighbor_lists = find_neighbors(faces, len(points))
    >>> label_pair_lists = sulcus_boundaries()
    >>> fraction_folds = 0.10  # low to speed up
    >>> min_sulcus_size = 50
    >>> sulcus_IDs, n_sulci = extract_sulci(label_pair_lists, labels,
    >>>     depth_file, area_file, neighbor_lists, fraction_folds,
    >>>     min_sulcus_size, do_fill_holes=False)
    >>> # Write results to vtk file and view with mayavi2:
    >>> rewrite_scalars(depth_file, 'test_extract_sulci.vtk',
    >>>                 sulcus_IDs, sulcus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')
    >>> # Write and view manual labels restricted to sulci:
    >>> rewrite_scalars(depth_file, 'test_extract_sulci_labels.vtk',
    >>>                 labels, sulcus_IDs)
    >>> os.system('mayavi2 -m Surface -d test_extract_sulci_labels.vtk &')

    """
    import numpy as np
    from time import time
    from mindboggle.extract.extract_folds import extract_all_folds
    from mindboggle.utils.mesh_operations import detect_boundaries,\
        segment, fill_holes

    #---------------------------------------------------------------------------
    # Load deep vertices
    #---------------------------------------------------------------------------
    folds, min_depth, n_deep_vertices = extract_all_folds(depth_file,
        area_file, fraction_folds)
    indices_folds = [i for i,x in enumerate(folds) if x == 1]
    if len(indices_folds):

        #---------------------------------------------------------------------------
        # Segment sulcus folds of a surface mesh
        #---------------------------------------------------------------------------
        print("  Segment surface mesh into sulcus folds...")
        t0 = time()

        # Array of sulcus IDs for fold vertices, initialized as -1.
        # Since we do not touch gyral vertices and vertices whose labels
        # are not in the label list, or vertices having only one label,
        # their sulcus IDs will remain -1.
        sulcus_IDs = -1 * np.ones(len(folds))

        # Find all label boundary pairs within the fold
        indices_boundaries, label_pairs, unique_label_pairs = detect_boundaries(
            indices_folds, labels, neighbor_lists)

        # Construct seed lists containing indices
        # to the vertices of each label boundary
        seed_lists = []
        for label_pair_list in label_pair_lists:
            seed_lists.append([x for i,x in enumerate(indices_boundaries)
                if np.sort(label_pairs[i]).tolist() in label_pair_list])

        # Segment into sets of vertices connected to label boundary seeds
        print("    Segment into separate label-pair regions...")
        t1 = time()
        sulcus_IDs = segment(indices_folds, neighbor_lists,
            seed_lists, min_sulcus_size, spread_same_labels=True,
            labels=labels, label_pair_lists=label_pair_lists)
        n_sulci = len([x for x in list(set(sulcus_IDs)) if x != -1])
        print("    ...Segmented {0} sulcus folds in {1:.2f} seconds".
              format(n_sulci, time() - t1))

    else:
        n_sulci = 0

    #---------------------------------------------------------------------------
    # Fill holes in folds
    #---------------------------------------------------------------------------
    if n_sulci > 0 and do_fill_holes:

        # Find fold vertices that have not yet been segmented
        # (because they weren't sufficiently deep)
        t2 = time()
        vertices_to_segment = [i for i,x in enumerate(sulcus_IDs) if x==-1]

        # Segment holes in the folds
        print("    Segment holes...")
        holes = segment(vertices_to_segment, neighbor_lists,
            seed_lists=[], min_region_size=1,
            spread_same_labels=False, labels=[], label_pair_lists=[])
        n_holes = len([x for x in list(set(holes)) if x != -1])

        # If there are any holes
        if n_holes > 0:

            # Ignore the largest hole (the background) and renumber holes
            max_hole_size = 0
            max_hole_index = 0
            for ihole in range(n_holes):
                I = np.where(holes == ihole)
                if len(I) > max_hole_size:
                    max_hole_size = len(I)
                    max_hole_index = ihole
            holes[holes == max_hole_index] = -1
            if max_hole_index < n_holes:
                holes[holes > max_hole_index] -= 1
            n_holes -= 1
            print('    ...{0} holes segmented ({1:.2f} seconds)'.
                  format(n_holes, time() - t2))

            # Fill holes
            t3 = time()
            sulcus_IDs = fill_holes(sulcus_IDs, holes, n_holes, neighbor_lists)
            print("    Filled holes ({0:.2f} seconds)".format(time() - t3))

    print("  ...Extracted and filled {0} sulcus folds in {1:.2f} seconds".
          format(n_sulci, time() - t0))

    return sulcus_IDs, n_sulci


# Example for extract_sulci
if __name__ == "__main__":

    import os
    from mindboggle.utils.io_vtk import load_scalar, rewrite_scalars
    from mindboggle.utils.mesh_operations import find_neighbors
    from mindboggle.extract.extract_folds import extract_sulci
    from mindboggle.info.sulcus_boundaries import sulcus_boundaries

    data_path = os.environ['MINDBOGGLE_DATA']
    label_file = os.path.join(data_path, 'subjects', 'MMRR-21-1',
                 'label', 'lh.labels.DKT25.manual.vtk')
    depth_file = os.path.join(data_path, 'measures',
                 '_hemi_lh_subject_MMRR-21-1', 'lh.pial.depth.vtk')
    area_file = os.path.join(data_path, 'measures',
                '_hemi_lh_subject_MMRR-21-1', 'lh.pial.area.vtk')

    points, faces, labels, n_vertices = load_scalar(label_file, True)
    neighbor_lists = find_neighbors(faces, len(points))
    label_pair_lists = sulcus_boundaries()

    fraction_folds = 0.10  # low to speed up
    min_sulcus_size = 50

    sulcus_IDs, n_sulci = extract_sulci(label_pair_lists, labels,
        depth_file, area_file, neighbor_lists, fraction_folds,
        min_sulcus_size, do_fill_holes=False)

    # Write results to vtk file and view with mayavi2:
    rewrite_scalars(depth_file, 'test_extract_sulci.vtk',
                    sulcus_IDs, sulcus_IDs)
    os.system('mayavi2 -m Surface -d test_extract_sulci.vtk &')
