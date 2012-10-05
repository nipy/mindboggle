#!/usr/bin/python
"""
Identify sulci from surface files with fold IDs and labels.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
#import sys
#import numpy as np
#from utils.io_vtk import load_scalar, write_scalars
#from info.sulcus_boundaries import sulcus_boundaries
#from utils.mesh_operations import segment, detect_boundaries, \
#    find_neighbors, compute_distance

verbose = 1


def identify(labels, folds, neighbor_lists, sulcus_names, label_pair_lists):
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

    """
    import numpy as np
    from utils.mesh_operations import segment, detect_boundaries, compute_distance

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
                subfolds, n_subfolds = segment(unassigned, neighbor_lists,
                    seed_lists, min_region_size=50,
                    spread_same_labels=False, labels=[], label_pair_lists=[])

                # Assign a sulcus ID to each segment
                print("           Assign sulcus IDs to segmented vertices")
                for ifold in range(n_subfolds):
                    subfold = [i for i,x in enumerate(subfolds) if x == ifold]
                    sulcus_IDs[subfold] = IDs_remainder_pairs[ifold]

    return sulcus_IDs


if __name__ == "__main__":

    """
    run extract/identify_sulci.py /Applications/freesurfer/subjects/MMRR-21-1/label/lh.labels.DKT25.manual.vtk /Users/arno/Desktop/output/results/features/_hemi_lh_subject_MMRR-21-1/sulci.vtk /Users/arno/Documents/Projects/Mindboggle/mindboggle/mindboggle/info/sulcus_names.txt test_identify.vtk
    """
    import sys
    from utils.io_vtk import load_scalar, rewrite_scalars
    from utils.io_file import read_columns
    from info.sulcus_boundaries import sulcus_boundaries
    from utils.mesh_operations import find_neighbors

    # Arguments
    labels_file = sys.argv[1]
    folds_file = sys.argv[2]
    sulcus_names_file = sys.argv[3]
    vtk_file = sys.argv[4]

    # Load labels and folds (the second surface has to be inflated).
    points, faces, labels, n_vertices = load_scalar(labels_file, return_arrays=0)
    points, faces, fold_IDs, n_vertices = load_scalar(folds_file, return_arrays=0)
    fid = open(sulcus_names_file, 'r')
    sulcus_names = fid.readlines()
    sulcus_names = [x.strip('\n') for x in sulcus_names]

    # Calculate neighbor lists for all vertices
    neighbor_lists = find_neighbors(faces, len(points))

    # Prepare list of all unique sorted label pairs in the labeling protocol
    label_pair_lists = sulcus_boundaries()

    # Create a list of lists of folds
    unique_fold_IDs = set(fold_IDs)
    folds = []
    for id in unique_fold_IDs:
        if id > 0:
            fold = [i for i,x in enumerate(fold_IDs) if x == id]
            folds.append(fold)

    # Identify sulci from folds
    sulcus_IDs = identify(labels, folds, neighbor_lists, sulcus_names,
                          label_pair_lists)

    # Finally, write points, faces and sulcus_IDs to a new vtk file
    rewrite_scalars(labels_file, vtk_file, sulcus_IDs, filter_scalars=sulcus_IDs)
